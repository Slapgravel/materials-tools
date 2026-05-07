# Crystal Packing Analysis — Project Spec

**Goal:** Mine the CSD (~1.24M structures) to discover how molecular features
(functional groups, shape, rotatable bonds, etc.) influence crystal packing
motifs. Build a queryable database that can answer questions like:
*"What is the role of halogens in high-symmetry molecule packing?"*

**User:** Jason Potticary (materials scientist, Meta MSI team)

---

## Scope & Filters

- **In scope (Phase 1):** Unitary small organics only
  - `is_organic = True`
  - No organometallics, no complexes, no solvates (single-component: no `,` in formula)
  - MW < 400 amu
  - No polymeric structures
  - No disorder
  - Must have 3D coordinates
  - **Estimated count: ~156,000 structures**

- **Future expansion:** Designed to be upgradable to the full CSD (complexes,
  solvates, higher MW, organometallics) by relaxing filter predicates.

---

## Architecture

```
CC_PackingAnalysis/
├── CLAUDE.md              ← this file (project context for resuming)
├── crystal_packing.db     ← SQLite database (main deliverable)
├── 01_extract.py          ← Phase 1: filter + extract metadata & descriptors
├── 02_packing.py          ← Phase 2: compute packing geometry (neighbor shells)
├── 03_cluster.py          ← Phase 3: unsupervised motif classification
├── 04_query.py            ← Phase 4: queryable interface
└── progress.json          ← checkpoint file for resumable extraction
```

### Database Schema

**`structures`** — one row per CSD entry
- refcode (PK), formula, smiles, molecular_weight, density
- spacegroup, crystal_system, z_value, z_prime
- a, b, c, alpha, beta, gamma, cell_volume
- packing_coefficient
- is_centrosymmetric, is_sohncke

**`molecular_descriptors`** — one row per refcode
- refcode (FK), n_atoms, n_heavy_atoms, n_carbon, n_nitrogen, n_oxygen, n_sulfur
- n_fluorine, n_chlorine, n_bromine, n_iodine, n_other_hetero
- n_rings, n_aromatic_rings, n_fused_ring_systems
- n_rotatable_bonds
- n_hbond_donors, n_hbond_acceptors
- planarity (eigenvalue ratio), aspect_ratio_1, aspect_ratio_2
- pi_system_fraction (aromatic C / total C)

**`packing_neighbors`** — 14 rows per refcode (nearest neighbor shell)
- refcode (FK), neighbor_rank (1-14)
- centroid_distance, interplanar_angle
- plane_offset (slip distance along stacking direction)

**`packing_motifs`** — one row per refcode (filled by Phase 3)
- refcode (FK), cluster_id, motif_label
- descriptor vector used for clustering

**`functional_groups`** — multiple rows per refcode
- refcode (FK), group_name, count
- Groups: OH, NH2, COOH, NO2, CN, halogen, ether, ester, amide,
  sulfone, etc. (SMARTS-based detection)

---

## CSD API Details

- **Python:** `C:\Users\potticary\CCDC\ccdc-software\csd-python-api\miniconda\python.exe`
- **API version:** 3.2.0
- **Database size:** 1,241,941 entries
- Key APIs:
  - `io.EntryReader('CSD')` — iterate/index entries
  - `entry.molecule` — molecular object (atoms, bonds, rings, SMILES)
  - `entry.crystal` — crystal object (cell, SG, packing)
  - `crystal.packing_shell(n)` — n nearest molecular neighbors
  - `descriptors.MolecularDescriptors.atom_plane(*atoms)` — best-fit plane
  - `descriptors.MolecularDescriptors.atom_centroid(*atoms)` — centroid
  - `descriptors.GeometricDescriptors.point_distance(p1, p2)` — distance
  - `plane.plane_angle(other_plane)` — interplanar angle

### Packing geometry extraction (validated)

For each structure, compute packing_shell(14) and for each neighbor:
1. Centroid-centroid distance
2. Interplanar angle (0° = cofacial/stacked, 90° = edge-on)
3. Plane offset (slip along stacking direction)

**Validated on known structures:**
- Naphthalene (NAPHTA11): herringbone — neighbors at ~52° edge-on ✓
- Anthracene (ANTCEN): herringbone — neighbors at ~51° edge-on ✓
- Pyrene (PYRENE03): slipped-stack — cofacial neighbor at 3.89 Å ✓

### Timing estimates
- Filter scan (metadata only): ~2ms/entry → ~0.8 hours for full CSD
- Packing geometry: ~50-100ms/entry → ~4-8 hours for ~156K structures
- Total extraction: ~5-9 hours (resumable via checkpoint)

---

## Packing Motif Classification (Phase 3)

**Approach: Unsupervised clustering — no predefined categories.**

Feature vector per structure (from neighbor shell):
- Distance to 1st, 2nd, 3rd nearest neighbor
- Interplanar angles to 1st, 2nd, 3rd nearest
- Ratio of cofacial (angle < 20°) vs edge-on (angle > 60°) neighbors
- Mean/std of neighbor distances
- Mean/std of interplanar angles

Clustering method: HDBSCAN (density-based, finds natural clusters, handles
noise). Followed by PCA visualization and manual labeling of discovered motifs.

**Expected motif types** (for validation, not as input):
- Herringbone (T-shaped, ~50-60° angles)
- Slipped-stack / brickwork (cofacial, small angles)
- Sandwich / cofacial (face-to-face, ~0° angles)
- Gamma (intermediate angles)
- Disordered / amorphous-like (no clear pattern)

---

## Query Interface (Phase 4)

SQLite + Python CLI. User types a natural-language question, which gets
translated to SQL queries and statistical analysis. Examples:

- "What are the role of halogens in high-symmetry molecule packing?"
- "Show me all herringbone-packing molecules with fluorine"
- "How does planarity correlate with packing coefficient?"
- "Which functional groups favor cofacial stacking?"

---

## Progress Tracking

Current status: **Phase 0 — Architecture & validation complete**

- [x] CSD API validated (v3.2.0, 1.24M entries)
- [x] Packing geometry extraction validated on 3 structures
- [x] Filter criteria defined (~156K structures estimated)
- [x] Database schema designed
- [x] Project folder created
- [ ] Phase 1: Extract metadata + molecular descriptors → SQLite
- [ ] Phase 2: Compute packing geometry (neighbor shells)
- [ ] Phase 3: Unsupervised clustering of packing motifs
- [ ] Phase 4: Query interface
- [ ] Phase 5 (optional): Expand to full CSD

---

## How to Resume

Point Claude at this folder:
```
C:\Users\potticary\Documents\PySandbox\materials-tools\crystallography\CC_PackingAnalysis\CLAUDE.md
```
Claude will read this file, check `progress.json` and `crystal_packing.db`
for current state, and continue from where it left off.
