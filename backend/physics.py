"""
QuantumHelix physics engine.

Models the proton in a DNA base-pair hydrogen bond as a 1D symmetric
double-well potential:  V(x) = V0 * ((x/a)^2 - 1)^2  (+ optional tilt)

- wkb():            WKB tunneling probability + attempt frequency + thermal hopping
- score_sequence(): per-base-pair mutagenic tautomer propensity along a sequence
- simulate():       time-dependent Schrodinger equation (split-step Fourier)

NOTE: PAIR_PARAMS are illustrative values chosen to be in a physically
plausible range. Replace them with DFT-derived barriers from the literature
(e.g. Slocombe et al., PCCP 2021; Winokan et al., Sci. Rep. 2023) before
drawing any scientific conclusion.
"""
import math
import numpy as np

HBAR = 1.054571817e-34   # J s
MP = 1.67262192e-27      # kg (proton mass)
EV = 1.602176634e-19     # J
KB = 1.380649e-23        # J/K
ANG = 1e-10              # m

# Barrier height V0 (eV), half-distance between wells a (Angstrom),
# tautomer asymmetry (eV, how much higher the "wrong" well sits)
PAIR_PARAMS = {
    "GC": {"V0": 0.36, "a": 0.55, "asym": 0.10},
    "AT": {"V0": 0.42, "a": 0.55, "asym": 0.20},
}
# Illustrative stacking modifier: each neighbouring pair shifts the barrier (eV)
NEIGHBOR_SHIFT = {"AT": -0.010, "GC": 0.004}
COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}


def pair_type(base: str) -> str:
    return "GC" if base in "GC" else "AT"


def clean_sequence(seq: str) -> str:
    """Strip FASTA headers, whitespace and non-ACGT characters."""
    lines = [l for l in seq.splitlines() if not l.startswith(">")]
    return "".join(c for c in "".join(lines).upper() if c in "ACGT")


def wkb(V0_eV: float, a_A: float, mass_ratio: float = 1.0, T_K: float = 310.0, n: int = 400) -> dict:
    """WKB tunneling from the zero-point level of one well through the central barrier."""
    m = mass_ratio * MP
    V0 = V0_eV * EV
    a = a_A * ANG
    omega = math.sqrt(8 * V0 / (m * a * a))      # harmonic frequency at the well bottom
    E0 = HBAR * omega / 2                          # zero-point energy
    nu = omega / (2 * math.pi)                     # attempt frequency

    T = 1.0
    if E0 < V0:
        xt = a * math.sqrt(1 - math.sqrt(E0 / V0))  # classical turning points +/- xt
        x = np.linspace(-xt, xt, n + 1)
        V = V0 * ((x / a) ** 2 - 1) ** 2
        integrand = np.sqrt(np.clip(2 * m * (V - E0), 0, None))
        I = float(np.sum((integrand[1:] + integrand[:-1]) / 2 * np.diff(x)))
        T = math.exp(-2 * I / HBAR)

    k_tun = nu * T
    k_th = nu * math.exp(-max(0.0, V0 - E0) / (KB * T_K))
    return {
        "E0_eV": E0 / EV,
        "attempt_hz": nu,
        "T": T,
        "k_tunnel": k_tun,
        "k_thermal": k_th,
        "tunnel_share": k_tun / (k_tun + k_th),
    }


def score_sequence(seq: str, T_K: float = 310.0) -> dict:
    seq = clean_sequence(seq)
    if not seq:
        raise ValueError("Sequence contains no A, C, G or T bases.")

    sites = []
    for i, b in enumerate(seq):
        pt = pair_type(b)
        p = PAIR_PARAMS[pt]
        V0 = p["V0"]
        for j in (i - 1, i + 1):
            if 0 <= j < len(seq):
                V0 += NEIGHBOR_SHIFT[pair_type(seq[j])]
        r = wkb(V0, p["a"], 1.0, T_K)
        f_taut = 1 / (1 + math.exp(p["asym"] * EV / (KB * T_K)))  # equilibrium tautomer fraction
        raw = (r["k_tunnel"] + r["k_thermal"]) * f_taut
        cpg = (b == "C" and i + 1 < len(seq) and seq[i + 1] == "G") or \
              (b == "G" and i > 0 and seq[i - 1] == "C")
        sites.append({
            "i": i, "base": b, "comp": COMP[b], "pair": pt,
            "V0": round(V0, 4), "logp": math.log10(raw),
            "share": r["tunnel_share"], "k_tunnel": r["k_tunnel"], "k_thermal": r["k_thermal"],
            "cpg": cpg,
        })

    lo = min(s["logp"] for s in sites)
    hi = max(s["logp"] for s in sites)
    for s in sites:
        s["score"] = (s["logp"] - lo) / (hi - lo) if hi > lo else 0.5

    gc = sum(1 for b in seq if b in "GC")
    top = sorted(sites, key=lambda s: -s["logp"])[:5]
    peak = top[0]
    return {
        "sequence": seq,
        "temperature": T_K,
        "sites": sites,
        "stats": {
            "length": len(seq),
            "gc_percent": 100 * gc / len(seq),
            "hotspots": sum(1 for s in sites if s["score"] >= 0.85),
            "mean_tunnel_share": sum(s["share"] for s in sites) / len(sites),
            "cpg_sites": sum(1 for s in sites if s["cpg"]),
            "peak_index": peak["i"],
        },
        "top": [s["i"] for s in top],
    }


def simulate(V0_eV=0.10, a_A=0.55, asym_eV=0.01, mass_ratio=1.0, total_fs=2000.0, samples=200) -> dict:
    """Split-step Fourier solution of the TDSE in units hbar=1, m=proton mass, length=Angstrom."""
    EU = HBAR ** 2 / (MP * ANG ** 2) / EV       # energy unit in eV (~4.15 meV)
    TU_fs = HBAR / (EU * EV) * 1e15             # time unit in fs (~158.7 fs)
    N, L = 512, 5.0
    x = np.linspace(-L / 2, L / 2, N, endpoint=False)
    dx = L / N
    k = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    m = mass_ratio
    V0, asym = V0_eV / EU, asym_eV / EU
    a = a_A
    V = V0 * ((x / a) ** 2 - 1) ** 2 + asym * (x + a) / (2 * a)

    dt = 0.0015
    expV = np.exp(-1j * V * dt / 2)
    expK = np.exp(-1j * (k ** 2 / (2 * m)) * dt)

    omega = math.sqrt(8 * V0 / (m * a * a))
    psi = np.exp(-m * omega * (x + a) ** 2 / 2).astype(complex)
    psi /= np.sqrt(np.sum(np.abs(psi) ** 2) * dx)

    total_steps = int(total_fs / TU_fs / dt)
    every = max(1, total_steps // samples)
    times, p_right = [], []
    for s in range(total_steps + 1):
        if s % every == 0:
            times.append(s * dt * TU_fs)
            p_right.append(float(np.sum(np.abs(psi[x > 0]) ** 2) * dx))
        psi = expV * np.fft.ifft(expK * np.fft.fft(expV * psi))
    return {"times_fs": times, "p_right": p_right, "E0_eV": omega / 2 * EU}
