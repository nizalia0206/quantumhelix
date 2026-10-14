# QuantumHelix

Simulates proton tunneling in DNA base pairs and maps where a sequence is most vulnerable to tautomer-driven mutation.

## Files

| File | What it does |
|---|---|
| `frontend/index.html` | The full app: 3D hero helix, live Schrödinger simulator, 3D sequence scanner. Works on its own. |
| `backend/physics.py` | WKB tunneling, sequence scoring, split-step TDSE solver (NumPy) |
| `backend/main.py` | FastAPI server: `/health`, `/basepair`, `/analyze`, `/analyze-file`, `/simulate` |
| `backend/requirements.txt` | Python dependencies |

## Run it on your laptop

```bash
# 1. Backend
cd backend
python -m venv venv
source venv/bin/activate          # Windows: venv\Scripts\activate
pip install -r requirements.txt
uvicorn main:app --reload --port 8000
# API docs: http://localhost:8000/docs

# 2. Frontend (new terminal)
cd frontend
python -m http.server 5500
# open http://localhost:5500
```

The badge in the top-right shows **Engine · Python API** when the page reaches the backend. If the backend is off, it says **Engine · Browser** and the same physics runs in JavaScript, so the demo never breaks.

To point the page at a deployed backend: `http://localhost:5500/?api=https://your-api.onrender.com`

## Deploy to Vercel (page + API on one URL)

The project root is already set up for Vercel:

- `index.html` is the page
- `api/index.py` + `api/_physics.py` is the FastAPI backend, served at `/api/...`
- `vercel.json` routes `/api/*` to the backend
- `requirements.txt` lists the Python packages

Steps: push the `quantumhelix` folder to GitHub, then on vercel.com click **Add New → Project**, import the repo and click **Deploy** (leave all settings default). The badge shows **Engine · Python API** when the backend responds. API docs live at `/api/docs`.

`backend/` is only for running locally with uvicorn.
