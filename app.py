from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Optional
import pandas as pd
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))
from Fusion_Cycle import model as FusionCycleModel

app = FastAPI(title="AgriSol Solubility API")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

fc_model = FusionCycleModel()

class PredictRequest(BaseModel):
    solute_smiles: str
    solvent_smiles: str
    temperature_k: float = 298.15
    solvent_density: Optional[float] = None

class ScreenRequest(BaseModel):
    solute_smiles: str
    solvents: List[dict]
    temperature_k: float = 298.15

@app.get("/")
def root():
    return {"status": "AgriSol API is running"}

@app.post("/predict")
def predict(req: PredictRequest):
    try:
        import numpy as np
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
                return {
                    "logS": None,
                    "solubility_mol_per_L": None,
                    "status": "out_of_domain",
                    "message": "Solute-solvent pair is outside model domain."
                }
            
            return {
                "logS": logS_value,
                "solubility_mol_per_L": float(10 ** logS_value),
                "status": "success"
            }
            
        except Exception as model_error:
            return {
                "logS": None,
                "solubility_mol_per_L": None,
                "status": "out_of_domain",
                "message": f"Model could not process this pair: {str(model_error)}"
            }
            
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/screen")
def screen(req: ScreenRequest):
    try:
        rows = []
        for solvent in req.solvents:
            rows.append({
                "solute_smiles_canonical": req.solute_smiles,
                "solvent_smiles_canonical": solvent["smiles"],
                "Temperature [K]": req.temperature_k,
                "solvent_density": solvent.get("density")
            })
        df = pd.DataFrame(rows)
        logS = fc_model.calculate_solubility(df)
        results = []
        for i, solvent in enumerate(req.solvents):
            results.append({
                "solvent_name": solvent["name"],
                "solvent_smiles": solvent["smiles"],
                "logS": float(logS.iloc[i]),
                "solubility_mol_per_L": float(10 ** logS.iloc[i])
            })
        results.sort(key=lambda x: x["logS"], reverse=True)
        return {"results": results, "status": "success"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
