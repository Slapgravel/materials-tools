r"""
Phase 1: Extract metadata + molecular descriptors from CSD into SQLite.
Resumable via checkpoint. Run with the CCDC Python:

    C:\Users\potticary\CCDC\ccdc-software\csd-python-api\miniconda\python.exe 01_extract.py
"""

import os
import sys
import json
import time
import sqlite3
import traceback
import numpy as np
from ccdc import io, descriptors

DB_PATH = os.path.join(os.path.dirname(__file__), "crystal_packing.db")
CHECKPOINT_PATH = os.path.join(os.path.dirname(__file__), "progress.json")
COMMIT_EVERY = 500
REPORT_EVERY = 1000

# ---------------------------------------------------------------------------
# Database setup
# ---------------------------------------------------------------------------

SCHEMA = """
CREATE TABLE IF NOT EXISTS structures (
    refcode TEXT PRIMARY KEY,
    formula TEXT,
    smiles TEXT,
    molecular_weight REAL,
    density REAL,
    spacegroup TEXT,
    spacegroup_number INTEGER,
    crystal_system TEXT,
    z_value REAL,
    z_prime REAL,
    a REAL, b REAL, c REAL,
    alpha REAL, beta REAL, gamma REAL,
    cell_volume REAL,
    packing_coefficient REAL,
    is_centrosymmetric INTEGER,
    is_sohncke INTEGER,
    color TEXT,
    habit TEXT,
    polymorph TEXT
);

CREATE TABLE IF NOT EXISTS molecular_descriptors (
    refcode TEXT PRIMARY KEY,
    n_atoms INTEGER,
    n_heavy_atoms INTEGER,
    n_carbon INTEGER,
    n_nitrogen INTEGER,
    n_oxygen INTEGER,
    n_sulfur INTEGER,
    n_fluorine INTEGER,
    n_chlorine INTEGER,
    n_bromine INTEGER,
    n_iodine INTEGER,
    n_other_hetero INTEGER,
    n_rings INTEGER,
    n_aromatic_rings INTEGER,
    n_fused_ring_systems INTEGER,
    n_rotatable_bonds INTEGER,
    n_hbond_donors INTEGER,
    n_hbond_acceptors INTEGER,
    planarity REAL,
    aspect_ratio_1 REAL,
    aspect_ratio_2 REAL,
    pi_system_fraction REAL,
    FOREIGN KEY (refcode) REFERENCES structures(refcode)
);

CREATE TABLE IF NOT EXISTS packing_neighbors (
    refcode TEXT,
    neighbor_rank INTEGER,
    centroid_distance REAL,
    interplanar_angle REAL,
    PRIMARY KEY (refcode, neighbor_rank),
    FOREIGN KEY (refcode) REFERENCES structures(refcode)
);

CREATE TABLE IF NOT EXISTS packing_motifs (
    refcode TEXT PRIMARY KEY,
    cluster_id INTEGER,
    motif_label TEXT,
    FOREIGN KEY (refcode) REFERENCES structures(refcode)
);

CREATE TABLE IF NOT EXISTS extraction_errors (
    refcode TEXT PRIMARY KEY,
    phase TEXT,
    error TEXT
);
"""


def init_db(path):
    conn = sqlite3.connect(path)
    conn.executescript(SCHEMA)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.commit()
    return conn


# ---------------------------------------------------------------------------
# Checkpoint
# ---------------------------------------------------------------------------

def load_checkpoint():
    if os.path.exists(CHECKPOINT_PATH):
        with open(CHECKPOINT_PATH) as f:
            return json.load(f)
    return {"phase1_index": 0, "phase1_done": False,
            "phase2_index": 0, "phase2_done": False}


def save_checkpoint(cp):
    with open(CHECKPOINT_PATH, "w") as f:
        json.dump(cp, f, indent=2)


# ---------------------------------------------------------------------------
# Filters
# ---------------------------------------------------------------------------

def passes_filters(entry):
    if not entry.is_organic:
        return False
    if entry.is_polymeric:
        return False
    if entry.has_disorder:
        return False
    if "," in entry.formula:
        return False
    mol = entry.molecule
    if mol is None:
        return False
    if mol.molecular_weight > 400:
        return False
    if not mol.all_atoms_have_sites:
        return False
    return True


# ---------------------------------------------------------------------------
# Descriptor extraction
# ---------------------------------------------------------------------------

def extract_structure(entry, crystal):
    sg_num = None
    try:
        sn = crystal.spacegroup_number_and_setting
        if sn:
            sg_num = sn[0]
    except Exception:
        pass

    return {
        "refcode": entry.identifier,
        "formula": entry.formula,
        "smiles": entry.molecule.smiles if entry.molecule else None,
        "molecular_weight": entry.molecule.molecular_weight if entry.molecule else None,
        "density": crystal.calculated_density,
        "spacegroup": crystal.spacegroup_symbol,
        "spacegroup_number": sg_num,
        "crystal_system": crystal.crystal_system if hasattr(crystal, 'crystal_system') else None,
        "z_value": crystal.z_value,
        "z_prime": crystal.z_prime,
        "a": crystal.cell_lengths.a,
        "b": crystal.cell_lengths.b,
        "c": crystal.cell_lengths.c,
        "alpha": crystal.cell_angles.alpha,
        "beta": crystal.cell_angles.beta,
        "gamma": crystal.cell_angles.gamma,
        "cell_volume": crystal.cell_volume,
        "packing_coefficient": crystal.packing_coefficient,
        "is_centrosymmetric": int(crystal.is_centrosymmetric) if hasattr(crystal, 'is_centrosymmetric') else None,
        "is_sohncke": int(crystal.is_sohncke) if hasattr(crystal, 'is_sohncke') else None,
        "color": entry.color,
        "habit": entry.habit if hasattr(entry, 'habit') else None,
        "polymorph": entry.polymorph if hasattr(entry, 'polymorph') else None,
    }


def extract_mol_descriptors(entry):
    mol = entry.molecule
    atoms = mol.atoms
    heavy = [a for a in atoms if a.atomic_symbol != "H"]

    elem_counts = {}
    for a in atoms:
        elem_counts[a.atomic_symbol] = elem_counts.get(a.atomic_symbol, 0) + 1

    n_carbon = elem_counts.get("C", 0)
    n_arom_c = sum(1 for a in atoms if a.atomic_symbol == "C"
                   and any(r.is_aromatic for r in mol.rings if a in r.atoms))

    # Planarity and aspect ratios from PCA of heavy atom coordinates
    planarity = 0.0
    ar1 = 1.0
    ar2 = 1.0
    coords_list = [(a.coordinates.x, a.coordinates.y, a.coordinates.z)
                   for a in heavy if a.coordinates]
    if len(coords_list) >= 3:
        coords = np.array(coords_list)
        centered = coords - coords.mean(axis=0)
        cov = np.cov(centered.T)
        eigvals = np.sort(np.linalg.eigvalsh(cov))[::-1]
        if eigvals[0] > 1e-10:
            planarity = 1.0 - eigvals[2] / eigvals[0]
            ar1 = eigvals[0] / eigvals[1] if eigvals[1] > 1e-10 else 99.0
            ar2 = eigvals[1] / eigvals[2] if eigvals[2] > 1e-10 else 99.0

    # Count fused ring systems
    n_fused = 0
    try:
        ring_atoms_sets = [set(r.atoms) for r in mol.rings]
        visited = [False] * len(ring_atoms_sets)
        for i in range(len(ring_atoms_sets)):
            if visited[i]:
                continue
            n_fused += 1
            stack = [i]
            while stack:
                cur = stack.pop()
                if visited[cur]:
                    continue
                visited[cur] = True
                for j in range(len(ring_atoms_sets)):
                    if not visited[j] and ring_atoms_sets[cur] & ring_atoms_sets[j]:
                        stack.append(j)
    except Exception:
        pass

    return {
        "refcode": entry.identifier,
        "n_atoms": len(atoms),
        "n_heavy_atoms": len(heavy),
        "n_carbon": n_carbon,
        "n_nitrogen": elem_counts.get("N", 0),
        "n_oxygen": elem_counts.get("O", 0),
        "n_sulfur": elem_counts.get("S", 0),
        "n_fluorine": elem_counts.get("F", 0),
        "n_chlorine": elem_counts.get("Cl", 0),
        "n_bromine": elem_counts.get("Br", 0),
        "n_iodine": elem_counts.get("I", 0),
        "n_other_hetero": sum(v for k, v in elem_counts.items()
                              if k not in ("C", "H", "N", "O", "S",
                                           "F", "Cl", "Br", "I")),
        "n_rings": len(mol.rings),
        "n_aromatic_rings": sum(1 for r in mol.rings if r.is_aromatic),
        "n_fused_ring_systems": n_fused,
        "n_rotatable_bonds": sum(1 for b in mol.bonds if b.is_rotatable),
        "n_hbond_donors": sum(1 for a in atoms if a.is_donor),
        "n_hbond_acceptors": sum(1 for a in atoms if a.is_acceptor),
        "planarity": round(planarity, 6),
        "aspect_ratio_1": round(ar1, 4),
        "aspect_ratio_2": round(ar2, 4),
        "pi_system_fraction": round(n_arom_c / n_carbon, 4) if n_carbon > 0 else 0.0,
    }


# ---------------------------------------------------------------------------
# Packing geometry (Phase 2 — included here for single-pass option)
# ---------------------------------------------------------------------------

MD = descriptors.MolecularDescriptors
GD = descriptors.GeometricDescriptors


def extract_packing(entry, crystal):
    mol = entry.molecule
    heavy = [a for a in mol.atoms if a.atomic_symbol != "H" and a.coordinates]
    if len(heavy) < 3:
        return []

    ref_plane = MD.atom_plane(*heavy)
    ref_centroid = MD.atom_centroid(*heavy)

    try:
        shell = crystal.packing_shell(14)
    except Exception:
        return []

    neighbors = []
    for comp in shell.components:
        comp_heavy = [a for a in comp.atoms
                      if a.atomic_symbol != "H" and a.coordinates]
        if len(comp_heavy) < 3:
            continue
        comp_centroid = MD.atom_centroid(*comp_heavy)
        dist = GD.point_distance(ref_centroid, comp_centroid)
        if dist < 0.5:  # skip self-overlap
            continue
        comp_plane = MD.atom_plane(*comp_heavy)
        angle = abs(ref_plane.plane_angle(comp_plane))
        if angle > 90:
            angle = 180 - angle
        neighbors.append((dist, angle))

    neighbors.sort(key=lambda x: x[0])

    rows = []
    for rank, (dist, angle) in enumerate(neighbors[:14], start=1):
        rows.append({
            "refcode": entry.identifier,
            "neighbor_rank": rank,
            "centroid_distance": round(dist, 4),
            "interplanar_angle": round(angle, 2),
        })
    return rows


# ---------------------------------------------------------------------------
# SQL insert helpers
# ---------------------------------------------------------------------------

def insert_structure(conn, data):
    cols = list(data.keys())
    placeholders = ", ".join("?" for _ in cols)
    sql = f"INSERT OR IGNORE INTO structures ({', '.join(cols)}) VALUES ({placeholders})"
    conn.execute(sql, [data[c] for c in cols])


def insert_mol_descriptors(conn, data):
    cols = list(data.keys())
    placeholders = ", ".join("?" for _ in cols)
    sql = f"INSERT OR IGNORE INTO molecular_descriptors ({', '.join(cols)}) VALUES ({placeholders})"
    conn.execute(sql, [data[c] for c in cols])


def insert_packing(conn, rows):
    for data in rows:
        cols = list(data.keys())
        placeholders = ", ".join("?" for _ in cols)
        sql = f"INSERT OR IGNORE INTO packing_neighbors ({', '.join(cols)}) VALUES ({placeholders})"
        conn.execute(sql, [data[c] for c in cols])


def insert_error(conn, refcode, phase, error):
    conn.execute(
        "INSERT OR REPLACE INTO extraction_errors (refcode, phase, error) VALUES (?, ?, ?)",
        (refcode, phase, error[:500]))


# ---------------------------------------------------------------------------
# Main extraction loop
# ---------------------------------------------------------------------------

def run_extraction(skip_packing=False):
    conn = init_db(DB_PATH)
    cp = load_checkpoint()
    reader = io.EntryReader("CSD")
    total = len(reader)

    start_idx = cp.get("phase1_index", 0)
    if cp.get("phase1_done"):
        print("Phase 1 already complete. Delete progress.json to re-run.")
        return

    print(f"CSD: {total} entries")
    print(f"Resuming from index {start_idx}")
    print(f"Packing geometry: {'SKIP' if skip_packing else 'ON'}")
    print(f"Database: {DB_PATH}")
    print()

    n_passed = 0
    n_skipped = 0
    n_errors = 0
    t_start = time.time()
    t_last_report = t_start

    for idx in range(start_idx, total):
        try:
            entry = reader[idx]
        except Exception:
            n_errors += 1
            continue

        try:
            if not passes_filters(entry):
                n_skipped += 1
                if (idx + 1) % REPORT_EVERY == 0:
                    _report(idx, total, n_passed, n_skipped, n_errors, t_start)
                continue
        except Exception:
            n_skipped += 1
            continue

        refcode = entry.identifier
        try:
            crystal = entry.crystal
            struct_data = extract_structure(entry, crystal)
            mol_data = extract_mol_descriptors(entry)
            insert_structure(conn, struct_data)
            insert_mol_descriptors(conn, mol_data)

            if not skip_packing:
                packing_rows = extract_packing(entry, crystal)
                insert_packing(conn, packing_rows)

            n_passed += 1

        except Exception as e:
            n_errors += 1
            try:
                insert_error(conn, refcode, "extract", traceback.format_exc())
            except Exception:
                pass

        if (idx + 1) % COMMIT_EVERY == 0:
            conn.commit()
            cp["phase1_index"] = idx + 1
            save_checkpoint(cp)

        if (idx + 1) % REPORT_EVERY == 0:
            _report(idx, total, n_passed, n_skipped, n_errors, t_start)

    conn.commit()
    cp["phase1_index"] = total
    cp["phase1_done"] = True
    save_checkpoint(cp)
    conn.close()

    elapsed = time.time() - t_start
    print(f"\nDone. {n_passed} structures extracted, {n_errors} errors, "
          f"{n_skipped} skipped in {elapsed/3600:.1f} hours.")


def _report(idx, total, n_passed, n_skipped, n_errors, t_start):
    elapsed = time.time() - t_start
    pct = (idx + 1) / total * 100
    rate = (idx + 1) / elapsed if elapsed > 0 else 0
    eta_h = (total - idx - 1) / rate / 3600 if rate > 0 else 0
    print(f"[{pct:5.1f}%] idx={idx+1}/{total}  "
          f"extracted={n_passed}  skipped={n_skipped}  errors={n_errors}  "
          f"rate={rate:.0f}/s  ETA={eta_h:.1f}h")


# ---------------------------------------------------------------------------

if __name__ == "__main__":
    skip_pack = "--no-packing" in sys.argv
    run_extraction(skip_packing=skip_pack)
