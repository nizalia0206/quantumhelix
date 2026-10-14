"""QuantumHelix API — run with:  uvicorn main:app --reload --port 8000"""
from io import StringIO

from fastapi import FastAPI, File, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field

import physics

app = FastAPI(title="QuantumHelix API", version="0.1.0")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])

MAX_BP = 5000


class AnalyzeRequest(BaseModel):
    sequence: str
    temperature: float = Field(310.0, ge=200, le=400)


class SimulateRequest(BaseModel):
    V0: float = Field(0.10, gt=0.01, le=1.0, description="Barrier height (eV)")
    a: float = Field(0.55, ge=0.2, le=1.2, description="Half-distance between wells (Angstrom)")
    asym: float = Field(0.01, ge=0, le=0.3, description="Tautomer asymmetry (eV)")
    mass: float = Field(1.0, ge=1.0, le=3.0, description="1 = proton, 2 = deuteron")
    total_fs: float = Field(2000.0, gt=0, le=20000)


@app.get("/health")
def health():
    return {"status": "ok", "engine": "python"}


@app.get("/basepair")
def basepair(pair: str = "GC", temperature: float = 310.0, mass: float = 1.0):
    pair = pair.upper()
    if pair not in physics.PAIR_PARAMS:
        raise HTTPException(400, "pair must be GC or AT")
    p = physics.PAIR_PARAMS[pair]
    return {"pair": pair, "params": p, **physics.wkb(p["V0"], p["a"], mass, temperature)}


@app.post("/analyze")
def analyze(req: AnalyzeRequest):
    seq = physics.clean_sequence(req.sequence)
    if len(seq) > MAX_BP:
        raise HTTPException(400, f"Sequence too long: {len(seq)} bp (max {MAX_BP}).")
    try:
        return physics.score_sequence(seq, req.temperature)
    except ValueError as e:
        raise HTTPException(400, str(e))


@app.post("/analyze-file")
async def analyze_file(file: UploadFile = File(...), temperature: float = 310.0):
    from Bio import SeqIO
    text = (await file.read()).decode("utf-8", errors="ignore")
    records = list(SeqIO.parse(StringIO(text), "fasta"))
    seq = str(records[0].seq) if records else text
    return analyze(AnalyzeRequest(sequence=seq, temperature=temperature))


@app.post("/simulate")
def simulate(req: SimulateRequest):
    return physics.simulate(req.V0, req.a, req.asym, req.mass, req.total_fs)
