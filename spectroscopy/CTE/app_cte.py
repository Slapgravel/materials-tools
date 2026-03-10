from flask import Flask, request, jsonify, render_template_string
import gemmi
import numpy as np
import re

app = Flask(__name__)

_NUM_WITH_UNCERT_RE = re.compile(
    r"""^\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*(?:\(\d+\))?\s*$"""
)

def strip_uncertainty(value: str) -> str:
    if value is None:
        return ""
    s = str(value).strip()
    if s in ("", ".", "?"):
        return ""
    m = _NUM_WITH_UNCERT_RE.match(s)
    return m.group(1) if m else s

def get_value(block, tag: str) -> str:
    try:
        v = str(block.find_value(tag)).strip()
        if v in (".", "?"):
            return ""
        return v
    except Exception:
        return ""

def to_float(s: str):
    s = strip_uncertainty(s)
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None

def compute_volume(a, b, c, alpha, beta, gamma):
    if any(v is None for v in (a, b, c, alpha, beta, gamma)):
        return None
    try:
        cell = gemmi.UnitCell(a, b, c, alpha, beta, gamma)
        return float(cell.volume)
    except Exception:
        return None

def parse_cif_bytes(filename: str, data: bytes):
    doc = gemmi.cif.read_string(data.decode("utf-8", errors="replace"))
    out = []
    for block in doc:
        T = to_float(get_value(block, "_cell_measurement_temperature"))
        a = to_float(get_value(block, "_cell_length_a"))
        b = to_float(get_value(block, "_cell_length_b"))
        c = to_float(get_value(block, "_cell_length_c"))
        alpha = to_float(get_value(block, "_cell_angle_alpha"))
        beta  = to_float(get_value(block, "_cell_angle_beta"))
        gamma = to_float(get_value(block, "_cell_angle_gamma"))
        V = compute_volume(a, b, c, alpha, beta, gamma)

        if T is None or (a is None and b is None and c is None and V is None):
            continue

        out.append({
            "file": filename,
            "data_block": str(block.name),
            "T": T,
            "a": a, "b": b, "c": c, "V": V,
        })
    return out

def anchored_slope(T, y, Tref, yref):
    """
    Fit y = yref + m*(T - Tref) (least squares, anchored at reference).
    Returns m, or None if degenerate.
    """
    T = np.asarray(T, dtype=float)
    y = np.asarray(y, dtype=float)
    dT = T - Tref
    dy = y - yref
    denom = np.sum(dT * dT)
    if denom == 0:
        return None
    m = float(np.sum(dT * dy) / denom)
    return m

def compute_cte_with_reference(records, key, tmin, tmax, Tref):
    """
    In window [tmin,tmax], compute:
      - choose yref from data point with temperature closest to Tref (within all records)
      - anchored slope m in window
      - alpha = m / yref
      - also provide per-point secant alphas relative to (Tref,yref) for the window
    """
    # Find reference yref from closest temperature point that has this key
    candidates = [(abs(r["T"] - Tref), r) for r in records if r.get("T") is not None and r.get(key) is not None]
    if not candidates:
        return None
    _, refrec = min(candidates, key=lambda x: x[0])
    Tref_used = float(refrec["T"])
    yref = float(refrec[key])
    if yref == 0:
        return None

    # Collect window points
    window = [r for r in records if r.get("T") is not None and tmin <= r["T"] <= tmax and r.get(key) is not None]
    if len(window) < 2:
        return {
            "Tref": Tref_used, "yref": yref,
            "n": len(window),
            "anchored_slope_dy_dT": None,
            "alpha_1_per_K": None,
            "note": "Need at least 2 points in window."
        }

    T = [float(r["T"]) for r in window]
    y = [float(r[key]) for r in window]

    m = anchored_slope(T, y, Tref_used, yref)
    if m is None:
        return None

    alpha = m / yref

    # Per-point secant alphas (optional but useful)
    secants = []
    for Ti, yi in zip(T, y):
        if Ti == Tref_used:
            continue
        secants.append((yi - yref) / (yref * (Ti - Tref_used)))

    sec_mean = float(np.mean(secants)) if secants else None
    sec_std = float(np.std(secants)) if secants else None

    return {
        "Tref": Tref_used,
        "yref": yref,
        "n": len(window),
        "anchored_slope_dy_dT": m,
        "alpha_1_per_K": alpha,
        "secant_alpha_mean_1_per_K": sec_mean,
        "secant_alpha_std_1_per_K": sec_std,
    }

@app.get("/")
def index():
    return render_template_string(INDEX_HTML)

@app.post("/api/parse")
def api_parse():
    files = request.files.getlist("files")
    records = []
    for f in files:
        records.extend(parse_cif_bytes(f.filename, f.read()))
    records.sort(key=lambda r: r["T"])
    temps = sorted({float(r["T"]) for r in records if r.get("T") is not None})
    return jsonify({"records": records, "temperatures": temps})

@app.post("/api/cte_ref")
def api_cte_ref():
    payload = request.get_json(force=True)
    records = payload.get("records", [])
    tmin = payload.get("tmin")
    tmax = payload.get("tmax")
    Tref = payload.get("Tref")

    if tmin is None or tmax is None or Tref is None:
        return jsonify({"error": "tmin, tmax, and Tref are required"}), 400

    result = {
        "tmin": tmin,
        "tmax": tmax,
        "Tref": Tref,
        "a": compute_cte_with_reference(records, "a", tmin, tmax, Tref),
        "b": compute_cte_with_reference(records, "b", tmin, tmax, Tref),
        "c": compute_cte_with_reference(records, "c", tmin, tmax, Tref),
        "V": compute_cte_with_reference(records, "V", tmin, tmax, Tref),
    }
    return jsonify(result)

INDEX_HTML = r"""
<!doctype html>
<html>
<head>
  <meta charset="utf-8" />
  <title>CIF CTE with reference temperature</title>
  <script src="https://cdn.plot.ly/plotly-2.30.0.min.js"></script>
  <style>
    body { font-family: Arial, sans-serif; margin: 18px; }
    #plot { width: 100%; height: 580px; }
    .row { margin: 10px 0; }
    input[type="number"] { width: 110px; }
    table { border-collapse: collapse; margin-top: 10px; }
    th, td { border: 1px solid #ccc; padding: 6px 8px; }
    th { background: #f6f6f6; }
  </style>
</head>
<body>
  <h2>Upload CIFs → drag-select window → choose Tref → compute CTE (anchored)</h2>

  <div class="row">
    <input id="files" type="file" multiple />
    <button id="btnUpload">Upload & Plot</button>
  </div>

  <div class="row">
    <label>Tmin: <input id="tmin" type="number" step="0.1"></label>
    <label>Tmax: <input id="tmax" type="number" step="0.1"></label>

    <label style="margin-left:12px;">Reference temperature (Tref):</label>
    <select id="tref"></select>

    <button id="btnCompute" style="margin-left:12px;">Compute CTE</button>
    <span style="margin-left:10px;color:#444;">Tip: drag-select points to set Tmin/Tmax.</span>
  </div>

  <div id="plot"></div>

  <h3>CTE results (anchored at Tref)</h3>
  <div id="out"></div>

<script>
let RECORDS = [];

const PARAM_STYLE = {
  a: { color: "red",   name: "a (Å)"  },
  b: { color: "blue",  name: "b (Å)"  },
  c: { color: "green", name: "c (Å)"  },
  V: { color: "black", name: "V (Å³)" },
};

function fmtAlpha(x){
  if (x === null || x === undefined) return "n/a";
  const micro = x * 1e6;
  return `${x.toExponential(6)} ( ${micro.toFixed(3)} ×10^-6 /K )`;
}
function fmtNum(x, digits=6){
  if (x === null || x === undefined) return "n/a";
  return Number(x).toFixed(digits);
}

function setDefaultWindow(){
  const Ts = RECORDS.map(r => r.T).filter(t => Number.isFinite(t));
  if (!Ts.length) return;
  document.getElementById("tmin").value = Math.min(...Ts);
  document.getElementById("tmax").value = Math.max(...Ts);
}

function populateTref(temps){
  const sel = document.getElementById("tref");
  sel.innerHTML = "";
  for (const t of temps){
    const opt = document.createElement("option");
    opt.value = t;
    opt.textContent = t;
    sel.appendChild(opt);
  }
  // default: pick the minimum temperature (common choice)
  if (temps.length) sel.value = temps[0];
}

function buildTracesGroupedLegend(){
  // One trace per parameter (legend: a,b,c,V)
  const traces = [];
  for (const key of ["a","b","c","V"]){
    const xs = [];
    const ys = [];
    for (const r of RECORDS){
      const v = r[key];
      if (Number.isFinite(r.T) && v !== null && v !== undefined){
        xs.push(r.T);
        ys.push(v);
      }
    }
    if (!xs.length) continue;

    const style = PARAM_STYLE[key];
    traces.push({
      x: xs,
      y: ys,
      mode: "markers",
      name: style.name,
      marker: { size: 7, color: style.color },
      // Separate volume axis; still on same plot but right axis.
      yaxis: (key === "V") ? "y2" : "y1",
      legendgroup: key,
    });
  }
  return traces;
}

function plotAll(){
  const traces = buildTracesGroupedLegend();

  const layout = {
    title: "Drag-select points to define [Tmin,Tmax]. Choose Tref from dropdown.",
    xaxis: { title: "Temperature (K)" },
    yaxis:  { title: "Lattice parameters (Å)" },
    yaxis2: { title: "Volume (Å³)", overlaying: "y", side: "right" },
    legend: { orientation: "h" },
    margin: { t: 60 },
    dragmode: "select"
  };

  Plotly.newPlot("plot", traces, layout, {responsive:true, displaylogo:false})
    .then(gd => {
      gd.on("plotly_selected", async (ev) => {
        if (!ev || !ev.points || ev.points.length === 0) return;
        const Ts = ev.points.map(p => p.x).filter(x => Number.isFinite(x));
        if (Ts.length < 2) return;
        document.getElementById("tmin").value = Math.min(...Ts);
        document.getElementById("tmax").value = Math.max(...Ts);
        await computeCTE();
      });
    });
}

function renderTable(data){
  const rows = [
    ["a", "Å"],
    ["b", "Å"],
    ["c", "Å"],
    ["V", "Å³"],
  ];

  let html = `<table>
    <tr>
      <th>Parameter</th>
      <th>Unit</th>
      <th>Tref used (K)</th>
      <th>Xref</th>
      <th># points in window</th>
      <th>Anchored slope dX/dT</th>
      <th>α = (1/Xref) dX/dT (1/K)</th>
      <th>Mean(secant α) (1/K)</th>
      <th>Std(secant α) (1/K)</th>
    </tr>`;

  for (const [key, unit] of rows){
    const r = data[key];
    if (!r){
      html += `<tr>
        <td>${key}</td><td>${unit}</td>
        <td colspan="7">n/a</td>
      </tr>`;
      continue;
    }
    html += `<tr>
      <td>${key}</td>
      <td>${unit}</td>
      <td>${fmtNum(r.Tref, 2)}</td>
      <td>${fmtNum(r.yref, 6)}</td>
      <td>${r.n ?? "n/a"}</td>
      <td>${r.anchored_slope_dy_dT === null ? "n/a" : Number(r.anchored_slope_dy_dT).toExponential(6)}</td>
      <td>${fmtAlpha(r.alpha_1_per_K)}</td>
      <td>${fmtAlpha(r.secant_alpha_mean_1_per_K)}</td>
      <td>${fmtAlpha(r.secant_alpha_std_1_per_K)}</td>
    </tr>`;
  }
  html += `</table>`;
  document.getElementById("out").innerHTML = html;
}

async function uploadAndPlot(){
  const filesEl = document.getElementById("files");
  if (!filesEl.files.length){
    alert("Choose one or more CIF files.");
    return;
  }
  const fd = new FormData();
  for (const f of filesEl.files) fd.append("files", f);

  const resp = await fetch("/api/parse", { method: "POST", body: fd });
  const data = await resp.json();

  RECORDS = data.records || [];
  if (!RECORDS.length){
    alert("No usable records found. Need _cell_measurement_temperature and cell parameters.");
    return;
  }

  populateTref(data.temperatures || []);
  setDefaultWindow();
  plotAll();

  document.getElementById("out").innerHTML =
    `<div>Loaded ${RECORDS.length} data-block records. Select a window and a reference temperature (Tref).</div>`;
}

async function computeCTE(){
  if (!RECORDS.length){
    alert("Upload CIF files first.");
    return;
  }
  const tmin = parseFloat(document.getElementById("tmin").value);
  const tmax = parseFloat(document.getElementById("tmax").value);
  const Tref = parseFloat(document.getElementById("tref").value);

  if (!Number.isFinite(tmin) || !Number.isFinite(tmax) || tmin > tmax){
    alert("Enter a valid Tmin/Tmax window.");
    return;
  }
  if (!Number.isFinite(Tref)){
    alert("Choose a valid Tref.");
    return;
  }

  const resp = await fetch("/api/cte_ref", {
    method: "POST",
    headers: {"Content-Type":"application/json"},
    body: JSON.stringify({records: RECORDS, tmin, tmax, Tref})
  });

  const data = await resp.json();
  if (data.error){
    document.getElementById("out").textContent = data.error;
    return;
  }
  renderTable(data);
}

document.getElementById("btnUpload").addEventListener("click", uploadAndPlot);
document.getElementById("btnCompute").addEventListener("click", computeCTE);
</script>
</body>
</html>
"""

if __name__ == "__main__":
    app.run(debug=True)