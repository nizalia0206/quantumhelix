"""Tests for the sensitivity experiment. Run:  python -m pytest tests -q   (or: python tests/test_experiment.py)"""
import json, os, sys
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "experiments")); sys.path.insert(0, os.path.join(ROOT, "api"))
import sensitivity as sx  # noqa: E402
import _physics as ph  # noqa: E402

def test_baseline_replica_matches_physics_engine():
    for s in sx.make_sequences():
        mine = sx.site_logp(s["seq"], sx.BASELINE)
        ref = np.array([x["logp"] for x in ph.score_sequence(s["seq"])["sites"]])
        assert np.allclose(mine, ref, atol=1e-3), s["id"]

def test_metrics_match_scipy():
    from scipy import stats
    rng = np.random.default_rng(1)
    for _ in range(30):
        x = rng.integers(0, 6, 40).astype(float); y = rng.integers(0, 4, 40).astype(float)
        assert abs(sx.kendall_b(x, y) - stats.kendalltau(x, y, variant="b").statistic) < 1e-12
        assert abs(sx.spearman(x, y) - stats.spearmanr(x, y).statistic) < 1e-12

def test_metric_edge_cases():
    x = np.arange(10.0)
    assert sx.kendall_b(x, x) == 1.0 and sx.kendall_b(x, -x) == -1.0
    assert np.isnan(sx.kendall_b(np.ones(10), x))            # constant ranking -> undefined
    assert sx.jaccard(frozenset({1, 2}), frozenset({2, 3})) == 1 / 3

def test_sequences_are_synthetic_fixed_and_valid():
    a, b = sx.make_sequences(), sx.make_sequences()
    assert a == b and len(a) == 4
    for s in a:
        assert set(s["seq"]) <= set("ACGT") and len(s["seq"]) == 60 and s["id"].startswith("SYN-")
    gc = [sum(c in "GC" for c in s["seq"]) / 60 for s in a]
    assert gc[0] < gc[1] < gc[2]

def test_experiment_is_deterministic():
    o1, _ = sx.run(n_lhs=40); o2, _ = sx.run(n_lhs=40)
    assert json.dumps(o1, sort_keys=True) == json.dumps(o2, sort_keys=True)

def test_all_scanned_scenarios_have_positive_barriers():
    """Worst-case barrier per pair type over each range box must stay positive (and above the zero-point energy)."""
    for rng_ in (sx.RANGES, sx.STRESS_RANGES):
        nb = [min(sx.NS[t] * k for k in rng_["nb_scale"]) for t in ("AT", "GC")]   # most negative shift per neighbour type
        worst_nb = 2 * min(nb)
        gc = sx.PP["GC"]["V0"] * rng_["V0_scale"][0] + worst_nb
        at = sx.PP["AT"]["V0"] * rng_["V0_scale"][0] + min(0, rng_["AT_offset"][0]) + worst_nb
        assert min(gc, at) > 0.07, (min(gc, at))

def test_published_findings_hold_in_results_file():
    o = json.load(open(os.path.join(ROOT, "results", "sensitivity_results.json")))
    core, stress = o["studies"]["core"], o["studies"]["stress"]
    for g, s in zip(core["global_"], o["sequences"]):
        assert s["top_class"] == "GC"
        assert g["frac_tau_ge_stable"] == 1.0 and g["frac_same_top_class"] == 1.0   # rank order invariant in plausible ranges
        assert g["jac"]["median"] < 0.9                                              # but hotspot calls are not
    assert core["s1_rank"][0] == "nb_scale" and stress["s1_rank"][0] == "nb_scale"
    # at nb_scale = 0 the ranking is exactly the A-T vs G-C split
    for k, s in enumerate(o["sequences"]):
        assert abs(o["oat"]["nb_scale"]["tau"][k][0] - s["tau_vs_gc_indicator"]) < 1e-3
    assert any(g["frac_tau_negative"] > 0 for g in stress["global_"])               # rankings can break under stress
    assert all(f["first_flip_value"] is not None and f["first_flip_value"] < -0.1 for f in o["breaking"]["AT_offset"]["flip"])
    assert all(i["tau"] == 1.0 for i in o["isotope_check"])
    assert core["frac_scenarios_with_zpe_above_barrier"] == 0 and stress["frac_scenarios_with_zpe_above_barrier"] == 0

if __name__ == "__main__":
    for n, f in list(globals().items()):
        if n.startswith("test_"):
            f(); print("ok", n)
