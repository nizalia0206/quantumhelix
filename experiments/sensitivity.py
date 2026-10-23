"""
QuantumHelix sensitivity analysis: "How robust are predicted proton-tunnelling
hotspots when the physical assumptions of the model change?"

Exploratory computational study. It re-uses api/_physics.py (wkb, PAIR_PARAMS,
NEIGHBOR_SHIFT, ...) and replicates score_sequence() with overridable
parameters. At baseline the replica matches physics.score_sequence() exactly
(checked in tests/test_experiment.py).

Run:  python experiments/sensitivity.py            (about 10-30 s)
Out:  results/*.json|csv|txt  and  public/investigation_data.js
Everything is seeded, so the numbers are reproducible.
"""
import argparse, csv, hashlib, json, math, os, random, sys
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "api"))
import _physics as ph  # noqa: E402

VERSION = "1.0"
SEED = 20261021          # one seed drives sequences, LHS and pair sampling
N_LHS = 1000
N_PAIRS = 3000
HOTSPOT_CUT = 0.85       # same rule the app uses for "hotspots"
STABLE_TAU = 0.8         # arbitrary, stated "stable" threshold for summaries
ROUND = 9                # log10(rate) rounded to 1e-9 before ranking (kills float noise)

# --------------------------------------------------------------------------
# Parameters. Baseline values are READ from api/_physics.py (illustrative).
# --------------------------------------------------------------------------
PP, NS = ph.PAIR_PARAMS, ph.NEIGHBOR_SHIFT
assert PP["GC"]["a"] == PP["AT"]["a"]
BASELINE = dict(V0_scale=1.0, AT_offset=0.0, a=PP["GC"]["a"], asym_scale=1.0,
                T=310.0, nb_scale=1.0, mass=1.0)
# Ranges: wide enough to test the model's own sensitivity (about +-20-50 % of the
# illustrative baselines). They are NOT claimed to be literature-derived ranges.
RANGES = {
    "V0_scale":   (0.80, 1.20),   # multiplies both barrier heights (0.36 / 0.42 eV)
    "AT_offset":  (-0.06, 0.06),  # eV added to the A-T barrier only (-0.06 => same barrier as G-C)
    "a":          (0.45, 0.65),   # half-distance between wells, angstrom (both pair types)
    "asym_scale": (0.50, 1.50),   # multiplies tautomer asymmetry energies (0.10 / 0.20 eV)
    "T":          (290.0, 330.0), # kelvin
    "nb_scale":   (0.0, 3.0),     # multiplies the neighbour (stacking) barrier shifts
}
# Stress ranges: deliberately wider, incl. reversed neighbour effect (nb_scale < 0) and a
# much lower A-T barrier. They probe where the ranking can break; some corners leave the
# regime where the barrier exceeds the zero-point energy (counted and reported).
STRESS_RANGES = {
    "V0_scale":   (0.70, 1.30),
    "AT_offset":  (-0.15, 0.06),
    "a":          (0.45, 0.65),
    "asym_scale": (0.25, 2.00),
    "T":          (280.0, 340.0),
    "nb_scale":   (-1.0, 3.0),
}
PARAMS = list(RANGES)
LABELS = {
    "V0_scale": "Barrier height (both pairs, x)",
    "AT_offset": "A-T barrier offset (eV)",
    "a": "Barrier half-width a (A)",
    "asym_scale": "Tautomer asymmetry (x)",
    "T": "Temperature (K)",
    "nb_scale": "Neighbour shift (x)",
}
OAT_VALUES = {
    "V0_scale": [0.8, 0.9, 1.0, 1.1, 1.2],
    "AT_offset": [-0.06, -0.03, 0.0, 0.03, 0.06],
    "a": [0.45, 0.50, 0.55, 0.60, 0.65],
    "asym_scale": [0.5, 0.75, 1.0, 1.25, 1.5],
    "T": [290.0, 300.0, 310.0, 320.0, 330.0],
    "nb_scale": [0.0, 0.5, 1.0, 2.0, 3.0],
}

# --------------------------------------------------------------------------
# Synthetic sequences (NOT real genes). Fixed seed, fixed length.
# --------------------------------------------------------------------------
def synth(n, gc, rng):
    out = []
    for _ in range(n):
        out.append(rng.choice("GC") if rng.random() < gc else rng.choice("AT"))
    return "".join(out)

def make_sequences(seed=SEED, n=60):
    rng = random.Random(seed)
    seqs = [
        ("SYN-AT25", "Synthetic, AT-rich (target 25% GC)", synth(n, 0.25, rng)),
        ("SYN-GC50", "Synthetic, balanced (target 50% GC)", synth(n, 0.50, rng)),
        ("SYN-GC75", "Synthetic, GC-rich (target 75% GC)", synth(n, 0.75, rng)),
        ("SYN-BLOCK", "Synthetic, 30 bp AT-rich block (20% GC) then 30 bp GC-rich block (80% GC)",
         synth(n // 2, 0.20, rng) + synth(n - n // 2, 0.80, rng)),
    ]
    return [dict(id=i, label=l, seq=s) for i, l, s in seqs]

# --------------------------------------------------------------------------
# Model (replica of physics.score_sequence with overridable parameters)
# --------------------------------------------------------------------------
_cache = {}
def _wkb(V0, a, m, T):
    key = (round(V0, 12), a, m, T)
    r = _cache.get(key)
    if r is None:
        r = _cache[key] = ph.wkb(V0, a, m, T)
    return r

def site_logp(seq, P, bad=None):
    """log10 of (k_tunnel + k_thermal) * tautomer fraction for every base pair.
    If `bad` (a list) is given, the number of sites whose zero-point energy is
    >= the barrier (V0 <= E0: model regime not meaningful) is appended to it."""
    out = np.empty(len(seq)); nbad = 0
    for i, b in enumerate(seq):
        pt = ph.pair_type(b)
        V0 = PP[pt]["V0"] * P["V0_scale"] + (P["AT_offset"] if pt == "AT" else 0.0)
        for j in (i - 1, i + 1):
            if 0 <= j < len(seq):
                V0 += NS[ph.pair_type(seq[j])] * P["nb_scale"]
        r = _wkb(V0, P["a"], P["mass"], P["T"])
        nbad += r["E0_eV"] >= V0
        f = 1 / (1 + math.exp(PP[pt]["asym"] * P["asym_scale"] * ph.EV / (ph.KB * P["T"])))
        out[i] = math.log10((r["k_tunnel"] + r["k_thermal"]) * f)
    if bad is not None:
        bad.append(int(nbad))
    return np.round(out, ROUND)

# --------------------------------------------------------------------------
# Ranking-stability metrics (tie-aware: the model produces very few distinct values)
# --------------------------------------------------------------------------
def avg_rank(x):
    _, inv, cnt = np.unique(x, return_inverse=True, return_counts=True)
    cum = np.cumsum(cnt)
    return (cum - (cnt - 1) / 2.0)[inv]

def spearman(x, y):
    rx, ry = avg_rank(x), avg_rank(y)
    if rx.std() == 0 or ry.std() == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])

def kendall_b(x, y):
    dx = np.sign(x[:, None] - x[None, :]); dy = np.sign(y[:, None] - y[None, :])
    nx, ny = (dx != 0).sum() / 2, (dy != 0).sum() / 2
    if nx == 0 or ny == 0:
        return float("nan")
    return float((dx * dy).sum() / 2 / math.sqrt(nx * ny))

def hotspot_set(lp):
    lo, hi = lp.min(), lp.max()
    if hi == lo:
        return frozenset()
    return frozenset(np.flatnonzero((lp - lo) / (hi - lo) >= HOTSPOT_CUT).tolist())

def jaccard(a, b):
    if not a and not b:
        return float("nan")
    return len(a & b) / len(a | b)

def top_class(lp, seq):
    gc = np.array([b in "GC" for b in seq])
    return "GC" if lp[gc].mean() > lp[~gc].mean() else ("AT" if lp[~gc].mean() > lp[gc].mean() else "tie")

def compare(base_lp, lp, seq):
    return dict(tau=kendall_b(base_lp, lp), rho=spearman(base_lp, lp),
                jac=jaccard(hotspot_set(base_lp), hotspot_set(lp)),
                same_top=top_class(base_lp, seq) == top_class(lp, seq),
                top=top_class(lp, seq), n_hot=len(hotspot_set(lp)))

# --------------------------------------------------------------------------
def lhs(n, rng, ranges):
    u = np.empty((n, len(PARAMS)))
    for k in range(len(PARAMS)):
        u[:, k] = (rng.permutation(n) + rng.random(n)) / n
    lo = np.array([ranges[p][0] for p in PARAMS]); hi = np.array([ranges[p][1] for p in PARAMS])
    return lo + u * (hi - lo)

def s1_binned(x, y, bins=10):
    """First-order sensitivity index estimate: Var(E[y|x]) / Var(y), x split in quantile bins.
    Returns 0 when y does not vary at all (nothing to explain)."""
    ok = ~np.isnan(y)
    x, y = x[ok], y[ok]
    if len(y) == 0 or y.var() == 0:
        return 0.0
    edges = np.quantile(x, np.linspace(0, 1, bins + 1)); idx = np.clip(np.searchsorted(edges, x, side="right") - 1, 0, bins - 1)
    used = [b for b in range(bins) if (idx == b).any()]
    means = np.array([y[idx == b].mean() for b in used]); w = np.array([(idx == b).sum() for b in used])
    return float(np.sum(w * (means - y.mean()) ** 2) / len(y) / y.var())

def r(x, d=4):
    return None if x is None or (isinstance(x, float) and math.isnan(x)) else round(float(x), d)

def quant(a):
    a = np.asarray(a, float); a = a[~np.isnan(a)]
    if len(a) == 0:
        return {}
    return dict(min=r(a.min()), p05=r(np.quantile(a, .05)), p25=r(np.quantile(a, .25)), median=r(np.median(a)),
                p75=r(np.quantile(a, .75)), p95=r(np.quantile(a, .95)), max=r(a.max()), mean=r(a.mean()))

def study(ranges, n_lhs, seed, seqs, base, base_lp):
    """Global Latin-hypercube study over `ranges`; every scenario is compared with the baseline."""
    rng = np.random.default_rng(seed)
    nseq = len(seqs)
    X = lhs(n_lhs, rng, ranges)
    res = {k: dict(tau=[], rho=[], jac=[], same=[], top=[]) for k in range(nseq)}
    LP = [np.empty((n_lhs, len(s["seq"]))) for s in seqs]
    hot = [np.zeros(len(s["seq"])) for s in seqs]
    invalid = np.zeros(n_lhs, bool)
    for ri, row in enumerate(X):
        P = dict(base); P.update({p: float(row[j]) for j, p in enumerate(PARAMS)})
        for k, s in enumerate(seqs):
            bad = []
            lp = site_logp(s["seq"], P, bad); LP[k][ri] = lp
            invalid[ri] |= bad[0] > 0
            c = compare(base_lp[k], lp, s["seq"])
            for key, v in (("tau", c["tau"]), ("rho", c["rho"]), ("jac", c["jac"]), ("same", c["same_top"]), ("top", c["top"])):
                res[k][key].append(v)
            for i in hotspot_set(lp):
                hot[k][i] += 1
    glob, s1t, s1j = [], {p: [] for p in PARAMS}, {p: [] for p in PARAMS}
    for k, s in enumerate(seqs):
        tau = np.array(res[k]["tau"], float); jac = np.array(res[k]["jac"], float)
        hist, _ = np.histogram(tau[~np.isnan(tau)], bins=20, range=(-1, 1))
        jh, _ = np.histogram(jac[~np.isnan(jac)], bins=20, range=(0, 1))
        for j, p in enumerate(PARAMS):
            s1t[p].append(s1_binned(X[:, j], tau)); s1j[p].append(s1_binned(X[:, j], jac))
        tops = res[k]["top"]
        glob.append(dict(id=s["id"], tau=quant(tau), rho=quant(res[k]["rho"]), jac=quant(jac),
                         tau_hist=hist.tolist(), jac_hist=jh.tolist(),
                         frac_tau_ge_stable=r(np.nanmean(tau >= STABLE_TAU)), frac_tau_negative=r(np.nanmean(tau < 0)),
                         frac_jac_ge_half=r(np.nanmean(jac >= 0.5)),
                         frac_same_top_class=r(np.mean(res[k]["same"])),
                         frac_top_GC=r(np.mean([t == "GC" for t in tops])), frac_top_AT=r(np.mean([t == "AT" for t in tops])),
                         persistence=[r(v / n_lhs) for v in hot[k]]))
    s1t_o = {p: dict(per_seq=[r(v) for v in s1t[p]], mean=r(np.mean(s1t[p]))) for p in PARAMS}
    s1j_o = {p: dict(per_seq=[r(v) for v in s1j[p]], mean=r(np.mean(s1j[p]))) for p in PARAMS}
    ranked = sorted(PARAMS, key=lambda p: -(s1t_o[p]["mean"] + s1j_o[p]["mean"]))
    # scenario-vs-scenario (independent of the baseline, which is itself illustrative)
    prng = np.random.default_rng(seed + 1)
    pt = [[] for _ in seqs]; pj = [[] for _ in seqs]
    for _ in range(N_PAIRS):
        i, j = prng.choice(n_lhs, 2, replace=False)
        for k in range(nseq):
            pt[k].append(kendall_b(LP[k][i], LP[k][j])); pj[k].append(jaccard(hotspot_set(LP[k][i]), hotspot_set(LP[k][j])))
    pairwise = [dict(id=s["id"], tau=quant(pt[k]), jac=quant(pj[k]), frac_tau_ge_stable=r(np.nanmean(np.array(pt[k]) >= STABLE_TAU)),
                     frac_tau_negative=r(np.nanmean(np.array(pt[k]) < 0))) for k, s in enumerate(seqs)]
    # two-parameter grid over the two most influential parameters
    p1, p2 = ranked[0], ranked[1]
    g1 = np.linspace(*ranges[p1], 9); g2 = np.linspace(*ranges[p2], 9)
    grid = dict(x=p1, y=p2, xvals=[r(v, 4) for v in g1], yvals=[r(v, 4) for v in g2],
                tau=[[[None] * 9 for _ in range(9)] for _ in seqs], jac=[[[None] * 9 for _ in range(9)] for _ in seqs])
    for ai, v1 in enumerate(g1):
        for bi, v2 in enumerate(g2):
            P = dict(base); P[p1] = float(v1); P[p2] = float(v2)
            for k, s in enumerate(seqs):
                lp = site_logp(s["seq"], P); c = compare(base_lp[k], lp, s["seq"])
                grid["tau"][k][bi][ai] = r(c["tau"]); grid["jac"][k][bi][ai] = r(c["jac"])
    return dict(ranges={p: list(v) for p, v in ranges.items()}, n=n_lhs, global_=glob, s1_tau=s1t_o, s1_jac=s1j_o, s1_rank=ranked,
                pairwise=pairwise, grid=grid, frac_scenarios_with_zpe_above_barrier=r(invalid.mean()))

def breaking_scans(seqs, base, base_lp):
    """Extended one-at-a-time scans: how far must one assumption move before the ranking changes?"""
    out = {}
    for p, vals in (("AT_offset", [round(v, 3) for v in np.arange(-0.25, 0.0601, 0.01)]),
                    ("nb_scale", [round(v, 3) for v in np.arange(-2.0, 3.0001, 0.25)])):
        rows = dict(values=vals, tau=[[] for _ in seqs], jac=[[] for _ in seqs], top=[[] for _ in seqs], invalid=[[] for _ in seqs])
        for v in vals:
            P = dict(base); P[p] = v
            for k, s in enumerate(seqs):
                bad = []; lp = site_logp(s["seq"], P, bad); c = compare(base_lp[k], lp, s["seq"])
                rows["tau"][k].append(r(c["tau"])); rows["jac"][k].append(r(c["jac"])); rows["top"][k].append(c["top"]); rows["invalid"][k].append(bad[0] > 0)
        # first value (moving away from baseline) at which the top class flips
        flips = []
        for k in range(len(seqs)):
            t = rows["top"][k]; bv = base[p]
            cand = [v for v, tc in zip(vals, t) if tc != t[vals.index(bv)] ] if bv in vals else []
            flips.append(dict(first_flip_value=(float(min(cand, key=lambda v: abs(v - bv))) if cand else None)))
        rows["flip"] = flips
        out[p] = rows
    return out

def run(n_lhs=N_LHS, seed=SEED):
    seqs = make_sequences(seed)
    base = dict(BASELINE)
    base_lp = [site_logp(s["seq"], base) for s in seqs]
    seq_out = []
    for s, lp in zip(seqs, base_lp):
        gc_ind = np.array([b in "GC" for b in s["seq"]], float); lo, hi = lp.min(), lp.max()
        seq_out.append(dict(id=s["id"], label=s["label"], seq=s["seq"], length=len(s["seq"]), gc_percent=r(100 * gc_ind.mean(), 1),
                            base_score=[r(v, 4) for v in (lp - lo) / (hi - lo)], base_logp=[r(v, 4) for v in lp],
                            distinct_levels=int(len(np.unique(lp))), n_hotspots=len(hotspot_set(lp)), top_class=top_class(lp, s["seq"]),
                            tau_vs_gc_indicator=r(kendall_b(lp, gc_ind))))
    oat, oat_rows = {}, []
    for p in PARAMS:
        oat[p] = {"values": OAT_VALUES[p], "tau": [[] for _ in seqs], "jac": [[] for _ in seqs], "same_top": [[] for _ in seqs]}
        for v in OAT_VALUES[p]:
            P = dict(base); P[p] = v
            for k, s in enumerate(seqs):
                c = compare(base_lp[k], site_logp(s["seq"], P), s["seq"])
                oat[p]["tau"][k].append(r(c["tau"])); oat[p]["jac"][k].append(r(c["jac"])); oat[p]["same_top"][k].append(bool(c["same_top"]))
                oat_rows.append([p, v, s["id"], r(c["tau"]), r(c["rho"]), r(c["jac"]), bool(c["same_top"]), c["top"], c["n_hot"]])
    iso = []
    for k, s in enumerate(seqs):
        P = dict(base); P["mass"] = 2.0
        c = compare(base_lp[k], site_logp(s["seq"], P), s["seq"])
        iso.append(dict(id=s["id"], tau=r(c["tau"]), jac=r(c["jac"]), same_top=bool(c["same_top"])))
    studies = {"core": study(RANGES, n_lhs, seed, seqs, base, base_lp),
               "stress": study(STRESS_RANGES, n_lhs, seed + 7, seqs, base, base_lp)}
    out = dict(meta=dict(version=VERSION, seed=seed, n_lhs=n_lhs, n_pairs=N_PAIRS, hotspot_cut=HOTSPOT_CUT, stable_tau=STABLE_TAU,
                         baseline=BASELINE, labels=LABELS, params=PARAMS, pair_params=PP, neighbor_shift=NS, rounding_decimals=ROUND,
                         parameters_are_illustrative=True),
               sequences=seq_out, oat=oat, isotope_check=iso, studies=studies, breaking=breaking_scans(seqs, base, base_lp))
    return out, oat_rows

def summary_text(o):
    L = ["QuantumHelix sensitivity analysis (exploratory computational study)",
         f"seed={o['meta']['seed']}  LHS scenarios per study={o['meta']['n_lhs']}  scenario pairs={o['meta']['n_pairs']}",
         "Parameters are illustrative; results show mathematical sensitivity of this model, not biological validity.", ""]
    for s in o["sequences"]:
        L.append(f"{s['id']}: {s['length']} bp, GC={s['gc_percent']}%, distinct score levels={s['distinct_levels']}, baseline top class={s['top_class']}, hotspots={s['n_hotspots']}, tau-b(model, GC-indicator)={s['tau_vs_gc_indicator']}")
    for name, st in o["studies"].items():
        L.append(f"\n=== {name.upper()} study  (scenarios with zero-point energy above a barrier: {st['frac_scenarios_with_zpe_above_barrier']}) ===")
        L.append("ranges: " + ", ".join(f"{p}={v}" for p, v in st["ranges"].items()))
        for s, g, pw in zip(o["sequences"], st["global_"], st["pairwise"]):
            L.append(f"{s['id']}: tau-b vs baseline median={g['tau']['median']} [p05={g['tau']['p05']}, p95={g['tau']['p95']}], share>={o['meta']['stable_tau']}={g['frac_tau_ge_stable']}, share<0={g['frac_tau_negative']}")
            L.append(f"          hotspot Jaccard median={g['jac']['median']} [p05={g['jac']['p05']}, p95={g['jac']['p95']}], share>=0.5={g['frac_jac_ge_half']}; same top class={g['frac_same_top_class']}")
            L.append(f"          scenario-vs-scenario: tau-b median={pw['tau']['median']}, share>={o['meta']['stable_tau']}={pw['frac_tau_ge_stable']}; Jaccard median={pw['jac']['median']}")
        L.append("first-order sensitivity (binned S1, mean over sequences)  tau-b | Jaccard:")
        for p in st["s1_rank"]:
            L.append(f"   {p:11s} {st['s1_tau'][p]['mean']} | {st['s1_jac'][p]['mean']}")
        L.append(f"   2-D grid axes: {st['grid']['x']} x {st['grid']['y']}")
    L.append("\nOne-at-a-time minimum tau-b over the core sweep (per sequence):")
    for p in o["meta"]["params"]:
        L.append(f"   {p:11s} {[min(t for t in tt if t is not None) for tt in o['oat'][p]['tau']]}   min Jaccard {[min(t for t in tt if t is not None) for tt in o['oat'][p]['jac']]}")
    for p, b in o["breaking"].items():
        L.append(f"\nBreaking scan {p}: first top-class flip away from baseline = {[f['first_flip_value'] for f in b['flip']]}; any zero-point>barrier in scan: {any(any(x) for x in b['invalid'])}")
    L.append(f"\nIsotope check (mass x2): {o['isotope_check']}")
    return "\n".join(L)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=N_LHS); ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--out", default=os.path.join(ROOT, "results"))
    a = ap.parse_args()
    o, rows = run(a.n, a.seed)
    os.makedirs(a.out, exist_ok=True)
    o["meta"]["sha256_of_results"] = hashlib.sha256(json.dumps(o, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    with open(os.path.join(a.out, "sensitivity_results.json"), "w") as f:
        json.dump(o, f, sort_keys=True, indent=1)
    with open(os.path.join(a.out, "oat_sweeps.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["parameter", "value", "sequence", "kendall_tau_b", "spearman", "hotspot_jaccard", "same_top_class", "top_class", "n_hotspots"]); w.writerows(rows)
    txt = summary_text(o) + f"\n\nsha256 of results: {o['meta']['sha256_of_results']}\n"
    with open(os.path.join(a.out, "summary.txt"), "w") as f:
        f.write(txt)
    if os.path.abspath(a.out) == os.path.join(ROOT, "results"):
        with open(os.path.join(ROOT, "public", "investigation_data.js"), "w") as f:
            f.write("/* generated by experiments/sensitivity.py - do not edit by hand */\nwindow.QH_INVESTIGATION=" + json.dumps(o, sort_keys=True, separators=(",", ":")) + ";\n")
    print(txt)

if __name__ == "__main__":
    main()
