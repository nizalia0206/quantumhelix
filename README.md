# QuantumHelix

Simulates proton tunnelling in DNA base pairs and ranks the sites of a sequence by a simplified, **illustrative** model. It is a teaching and exploration tool, not a validated mutation predictor.

## Structure

```
public/index.html              the app (3D helix, live Schrödinger simulator, scanner, investigation section)
public/investigation_data.js   generated results the investigation section draws (do not edit by hand)
api/index.py                   FastAPI backend, routes under /api
api/_physics.py                WKB tunnelling, sequence scoring, TDSE solver (unchanged by the investigation)
experiments/sensitivity.py     the reproducible sensitivity analysis
tests/test_experiment.py       7 checks on the experiment and its published findings
results/                       generated: sensitivity_results.json, oat_sweeps.csv, summary.txt
requirements.txt, vercel.json  Python packages, Vercel config
```

## Run locally

Backend (optional, the page falls back to the same physics in JavaScript):
```
pip install -r requirements.txt
uvicorn api.index:app --reload --port 8000     # docs: http://localhost:8000/api/docs
```
Frontend:
```
python -m http.server 5500 --directory public   # open http://localhost:5500
```
The badge shows **Engine · Python API** when the backend is reachable, otherwise **Engine · Browser**. The 3D scene loads three.js from a CDN, so it needs internet unless you host that file yourself.

## Science Fair Mode

Click **Science Fair Mode** in the top bar, or open `index.html#fair` directly (handy for a booth bookmark). A full-screen, 5-step presentation (the question, three steps, what our test found) loops by itself with large text and a progress bar. It never navigates away from the page.

- Controls: Back, Pause/Play, Next, Restart, step dots, **Explore the experiment** (closes the overlay and jumps to the simulator) and **Exit** (Esc).
- Keys: Esc exit, ← → step, P or Space pause (Space on a focused button activates that button), R restart.
- Nothing essential depends on hover. Focus stays inside the overlay and returns to the button when it closes. The rest of the page is inert while it is open.
- With "reduce motion" set in the operating system, it starts paused with static drawings and does not auto-advance.
- The drawings are schematic: the double-well probability animation is illustrative, not computed from the simulator. A persistent strip on every step says so. The last step reads its numbers from `public/investigation_data.js`.

## The investigation: "Do the Hotspots Hold Up?"

Question: how robust are the predicted hotspots when the model's physical assumptions change? An exploratory computational study, not a pre-registered test. The hotspot-overlap metric and the "stress" study were added after the first run showed the rank order never changed.

Reproduce (about 10 seconds, fully seeded):
```
python experiments/sensitivity.py     # regenerates results/ and public/investigation_data.js
python tests/test_experiment.py       # runs the 7 checks
```
`results/summary.txt` ends with a SHA-256 of the results; two runs give the same hash.

Design: four **synthetic** 60 bp sequences (not real genes); one-at-a-time sweeps of six parameters; two 1,000-scenario Latin-hypercube studies (plausible ranges and wider stress ranges); 3,000 scenario-vs-scenario comparisons; Kendall τ-b, Spearman ρ, hotspot-set Jaccard overlap, top-class flips and a binned first-order sensitivity index. Scores are rounded to 1e-9 before ranking. The experiment re-implements `score_sequence` with overridable parameters; a test confirms it matches the physics engine at baseline.

Results actually obtained (seed 20261021):
- **Rank order is invariant in the plausible ranges.** τ-b ≥ 0.8 in 100% of scenarios for all four sequences. This follows from the model's built-in signs (G–C outranks A–T, A–T neighbours lower the barrier, G–C neighbours raise it), not from a discovery. Each sequence has only 6–8 distinct score levels; τ-b against a plain G/C indicator is 0.69–0.79.
- **The hotspot list is not stable.** Median overlap with the baseline hotspot set is 0.33–0.86 depending on the sequence. The app rescales scores to 0–1 and cuts at 0.85, so the set depends on score spread.
- **The neighbour-shift assumption matters most** (S1 = 0.62 for hotspot overlap, plausible study). Those shifts (−0.010 / +0.004 eV) have no literature source that we found. At zero shift the ranking is just A–T vs G–C.
- **Under stress the ranking can break:** up to 25% of scenarios have τ-b < 0 for the worst-affected sequence. A–T only outranks G–C when its barrier is lowered by about 0.18–0.20 eV.
- Doubling the mass (deuteron check) changes rates but not the ranking.

What this does **not** show: biological validity, mutation rates, disease risk, behaviour on real genomes, or any novel scientific result. Stability of a model's output under its own parameters is mathematical sensitivity, not evidence that the model is right.

## Limitations of the model

- One proton in a simplified 1-D double well; WKB tunnelling plus thermal hopping; no cellular environment, no repair or replication machinery.
- **The barrier, asymmetry and neighbour-shift values are illustrative**, not DFT- or experiment-derived. Parameter ranges in the study are ±20–50% of those baselines, not literature ranges.
- Site scores are relative rankings, not probabilities. "Hotspot" means highest relative score in this simulation.

## Comparison with published work (only where comparable)

- Slocombe, Al-Khalili & Sacchi (2021), *Phys. Chem. Chem. Phys.* 23, 4141–4150, https://doi.org/10.1039/D0CP05781A: DFT and machine-learning nudged-elastic-band study of A–T and G–C tautomerisation. A talk summary by one of the authors says lasting changes look unlikely for A–T but more significant for G–C. Our baseline also ranks G–C above A–T, but that order comes directly from the illustrative asymmetry values, so it is not independent validation.
- Slocombe, Sacchi & Al-Khalili (2022), *Commun. Phys.* 5, 109, https://doi.org/10.1038/s42005-022-00881-8: open-quantum-systems model of G–C; reports tunnelling outweighing the classical contribution. Our toy model gives about 5% tunnelling share for a 0.36 eV barrier. The physics differs, so the two are not comparable, and the mismatch is another reason to treat absolute rates as illustrative.
- Pushkaran & Arabi (2024), *Int. J. Biol. Macromol.* 277, 134051, https://doi.org/10.1016/j.ijbiomac.2024.134051: listed as a review for context; not used for any comparison here.

## Deploy (Vercel)

`vercel.json` rewrites `/api/*` to `api/index.py`; `public/` is served statically. After re-running the experiment, commit `public/investigation_data.js` and `results/`.
