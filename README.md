# QuantumHelix

Simulates proton tunneling in DNA base pairs and maps where a sequence is most vulnerable to tautomer-driven mutation.

## Structure

```
public/index.html     the app (3D helix, live Schrödinger simulator, 3D scanner)
api/index.py          FastAPI backend, routes under /api
api/_physics.py       WKB tunneling, sequence scoring, TDSE solver
requirements.txt      Python packages
vercel.json           Vercel config
```

## Run locally (Windows, from the project folder)

Terminal 1, backend:
```powershell
venv\Scripts\Activate.ps1
pip install -r requirements.txt
uvicorn api.index:app --reload --port 8000
```
API docs: http://localhost:8000/api/docs

Terminal 2, frontend:
```powershell
py -m http.server 5500 --directory public
```
Open http://localhost:5500. The badge shows **Engine · Python API** when connected, otherwise **Engine · Browser** (same physics in JavaScript).

## Deploy

Push to GitHub, import on vercel.com, Application Preset **Other**, Deploy. The page is served from `public/`, the API from `/api`.

Barrier values in `PAIR_PARAMS` are illustrative; replace with DFT values from the literature before drawing scientific conclusions.
