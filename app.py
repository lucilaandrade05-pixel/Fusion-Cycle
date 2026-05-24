from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Optional
import pandas as pd
import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))

app = FastAPI(title="AgriSol Solubility API")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# Load experimental database from CSV
_db_path = os.path.join(os.path.dirname(__file__), 'experimental_db.csv')
EXPERIMENTAL_DB = pd.read_csv(_db_path) if os.path.exists(_db_path) else pd.DataFrame()
# Try to load Fusion-Cycle model
fc_model = None
try:
    from Fusion_Cycle import model as FusionCycleModel
    fc_model = FusionCycleModel()
    print("Fusion-Cycle model loaded successfully")
except Exception as e:
    print(f"Fusion-Cycle model not available: {e}")

class PredictRequest(BaseModel):
    solute_smiles: str
    solvent_smiles: str
    temperature_k: float = 298.15
    solvent_density: Optional[float] = None

class ScreenRequest(BaseModel):
    solute_smiles: str
    solvents: List[dict]
    temperature_k: float = 298.15

def lookup_experimental(solute_smiles, solvent_smiles, temperature_k):
    match = EXPERIMENTAL_DB[
        (EXPERIMENTAL_DB['solute_smiles'] == solute_smiles.strip()) &
        (EXPERIMENTAL_DB['solvent_smiles'] == solvent_smiles.strip())
    ]
    if len(match) > 0:
        row = match.iloc[0]
        return {
            "logS": float(row['logS']),
            "solubility_mol_per_L": float(10 ** row['logS']),
            "solubility_g_L": float(row['solubility_g_L']),
            "status": "success",
            "source": "experimental",
            "solute_name": str(row['solute_name']),
            "solvent_name": str(row['solvent_name'])
        }
    return None

@app.get("/")
def root():
    return {"status": "AgriSol API is running", "experimental_pairs": len(EXPERIMENTAL_DB), "model_available": fc_model is not None}

@app.post("/predict")
def predict(req: PredictRequest):
    exp = lookup_experimental(req.solute_smiles, req.solvent_smiles, req.temperature_k)
    if exp:
        return exp
    if fc_model is None:
        return {"logS": None, "solubility_mol_per_L": None, "status": "out_of_domain", "message": "No experimental data found and model unavailable for this pair."}
    try:
        df = pd.DataFrame([{
            "solute_smiles_canonical": req.solute_smiles,
            "solvent_smiles_canonical": req.solvent_smiles,
            "Temperature [K]": req.temperature_k,
            "solvent_density": req.solvent_density
        }])
        try:
            logS = fc_model.calculate_solubility(df)
            logS_value = float(logS.iloc[0])
            if np.isnan(logS_value) or np.isinf(logS_value):
                return {"logS": None, "solubility_mol_per_L": None, "status": "out_of_domain", "message": "Pair outside model domain.", "source": "model"}
            return {"logS": logS_value, "solubility_mol_per_L": float(10 ** logS_value), "status": "success", "source": "model"}
        except Exception as model_error:
            return {"logS": None, "solubility_mol_per_L": None, "status": "out_of_domain", "message": str(model_error), "source": "model"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/screen")
def screen(req: ScreenRequest):
    results = []
    for solvent in req.solvents:
        exp = lookup_experimental(req.solute_smiles, solvent["smiles"], req.temperature_k)
        if exp:
            results.append({
                "solvent_name": solvent["name"],
                "solvent_smiles": solvent["smiles"],
                "logS": exp["logS"],
                "solubility_mol_per_L": exp["solubility_mol_per_L"],
                "solubility_g_L": exp.get("solubility_g_L"),
                "status": "success",
                "source": "experimental"
            })
        elif fc_model is not None:
            try:
                df = pd.DataFrame([{
                    "solute_smiles_canonical": req.solute_smiles,
                    "solvent_smiles_canonical": solvent["smiles"],
                    "Temperature [K]": req.temperature_k,
                    "solvent_density": solvent.get("density")
                }])
                logS = fc_model.calculate_solubility(df)
                logS_value = float(logS.iloc[0])
                results.append({
                    "solvent_name": solvent["name"],
                    "solvent_smiles": solvent["smiles"],
                    "logS": logS_value,
                    "solubility_mol_per_L": float(10 ** logS_value),
                    "status": "success",
                    "source": "model"
                })
            except Exception:
                results.append({"solvent_name": solvent["name"], "solvent_smiles": solvent["smiles"], "logS": None, "status": "out_of_domain", "source": "model"})
        else:
            results.append({"solvent_name": solvent["name"], "solvent_smiles": solvent["smiles"], "logS": None, "status": "out_of_domain", "source": "none"})
    results.sort(key=lambda x: x["logS"] if x["logS"] is not None else -999, reverse=True)
    return {"results": results, "status": "success"}

