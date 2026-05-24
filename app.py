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

# Load experimental database
EXPERIMENTAL_DB = pd.DataFrame([
    Gerando código Python para o banco de dados completo...

EXPERIMENTAL_DB = pd.DataFrame([
    {"solute_smiles":"O=C1N(CC2=CC=CC=C2Cl)OCC1(C)C","solvent_smiles":"O=C(OCCOC=1C=CC=CC1)C","logS":0.66,"solubility_g_L":1100.0,"solute_name":"Clomazone","solvent_name":"TEX MEX","temperature_K":323.15},
    {"solute_smiles":"O=C1N(CC2=CC=CC=C2Cl)OCC1(C)C","solvent_smiles":"O=C1CCC2OC1OC2","logS":0.52,"solubility_g_L":800.0,"solute_name":"Clomazone","solvent_name":"Cireno (Cyrene)","temperature_K":298.15},
    {"solute_smiles":"O=C1N(N)C(SC)=NN=C1C(C)(C)C","solvent_smiles":"O=C(OCC)C(=O)C","logS":0.51,"solubility_g_L":700.0,"solute_name":"Metribuzin","solvent_name":"Ethyl pyruvate","temperature_K":323.15},
    {"solute_smiles":"OC(C1(CC1)Cl)(CC2=C(Cl)C=CC=C2)CN3N=CNC3=S","solvent_smiles":"CN1C(CCC1)=O","logS":0.5,"solubility_g_L":1100.0,"solute_name":"Prothioconazole","solvent_name":"NMP","temperature_K":323.15},
    {"solute_smiles":"O=C1N(N)C(SC)=NN=C1C(C)(C)C","solvent_smiles":"O=C1CCC2OC1OC2","logS":0.45,"solubility_g_L":600.0,"solute_name":"Metribuzin","solvent_name":"Cireno (Cyrene)","temperature_K":323.15},
    {"solute_smiles":"O=C1N(N)C(SC)=NN=C1C(C)(C)C","solvent_smiles":"N#CCC(=O)OCC","logS":0.45,"solubility_g_L":600.0,"solute_name":"Metribuzin","solvent_name":"Ethyl cyanoacetate","temperature_K":323.15},
    {"solute_smiles":"O=C1N(CC2=CC=CC=C2Cl)OCC1(C)C","solvent_smiles":"O=C(N(C)C)CCCCCCCCC","logS":0.42,"solubility_g_L":635.0,"solute_name":"Clomazone","solvent_name":"Surfom CE 8434","temperature_K":298.15},
    {"solute_smiles":"O=C1N(CC2=CC=CC=C2Cl)OCC1(C)C","solvent_smiles":"O=C(OCC)C=1C=CC=CC1","logS":0.41,"solubility_g_L":615.0,"solute_name":"Clomazone","solvent_name":"Ethyl Benzoate","temperature_K":298.15},
    {"solute_smiles":"O=C1N(CC2=CC=CC=C2Cl)OCC1(C)C","solvent_smiles":"CC(=O)OCC1=CC=CC=C1","logS":0.409207082,"solubility_g_L":615.0,"solute_name":"Clomazone","solvent_name":"Surfonic AG 1705","temperature_K":298.15},
    {"solute_smiles":"O=C1N(CC2=CC=CC=C2Cl)OCC1(C)C","solvent_smiles":"O=C(OCCCC)CCC(=O)C","logS":0.4,"solubility_g_L":605.0,"solute_name":"Clomazone","solvent_name":"Butyl Levulinate","temperature_K":298.15},
    {"solute_smiles":"O=C1N(CC2=CC=CC=C2Cl)OCC1(C)C","solvent_smiles":"OCC1OC(OC1)(C)CC","logS":0.4,"solubility_g_L":600.0,"solute_name":"Clomazone","solvent_name":"GBK","temperature_K":298.15},
    {"solute_smiles":"O=C1N(CC2=CC=CC=C2Cl)OCC1(C)C","solvent_smiles":"CN1C(CCC1)=O","logS":0.4,"solubility_g_L":600.0,"solute_name":"Clomazone","solvent_name":"NMP","temperature_K":298.15},
    {"solute_smiles":"CNC(=O)C(=NOC)C1=CC=CC=C1OC2=CC=CC=C2","solvent_smiles":"O=C(OCC)C(=O)C","logS":0.39,"solubility_g_L":700.0,"solute_name":"Metominostrobin","solvent_name":"Ethyl pyruvate","temperature_K":323.15},
    {"solute_smiles":"CC(C1CC1)C(CN2C=NC=N2)(C3=CC=C(C=C3)Cl)O","solvent_smiles":"O=C(OCC)C(=O)C","logS":0.38,"solubility_g_L":700.0,"solute_name":"Cyproconazole","solvent_name":"Ethyl pyruvate","temperature_K":323.15},
    {"solute_smiles":"O=[N+](/N=C1NCCN\1CC2=CC=C(Cl)N=C2)[O-]","solvent_smiles":"CN1C(CCC1)=O","logS":0.37,"solubility_g_L":600.0,"solute_name":"Imidacloprid","solvent_name":"NMP","temperature_K":298.15},
    {"solute_smiles":"O=[N+](/N=C1NCCN\1CC2=CC=C(Cl)N=C2)[O-]","solvent_smiles":"CN1C(CCC1)=O","logS":0.37,"solubility_g_L":600.0,"solute_name":"Imidacloprid","solvent_name":"NMP","temperature_K":323.15},
    {"solute_smiles":"O=C(OC)/C(C1=CC=CC=C1OC2=NC=NC(OC3=CC=CC=C3C#N)=C2)=C/OC","solvent_smiles":"CN1C(CCC1)=O","logS":0.35,"solubility_g_L":900.0,"solute_name":"Azoxystrobin","solvent_name":"NMP","temperature_K":323.15},
    {"solute_smiles":"OC(C(C)(C)C)(CCC1=CC=C(Cl)C=C1)CN2C=NC=N2","solvent_smiles":"CN1C(CCC1)=O","logS":0.324644741,"solubility_g_L":650.0,"solute_name":"Tebuconazole","solvent_name":"NMP","temperature_K":298.15},
    {"solute_smiles":"OC(C(C)(C)C)(CCC1=CC=C(Cl)C=C1)CN2C=NC=N2","solvent_smiles":"O=C(OCC)C(=O)C","logS":0.289882635,"solubility_g_L":600.0,"solute_name":"Tebuconazole","solvent_name":"Ethyl pyruvate","temperature_K":298.15},
    {"solute_smiles":"OC(C(C)(C)C)(CCC1=CC=C(Cl)C=C1)CN2C=NC=N2","solvent_smiles":"CCCCCCc1ccccc1","logS":0.289882635,"solubility_g_L":600.0,"solute_name":"Tebuconazole","solvent_name":"Solvesso 150","temperature_K":323.15},
    {"solute_smiles":"O=C(OC)/C(C1=CC=CC=C1OC2=NC=NC(OC3=CC=CC=C3C#N)=C2)=C/OC","solvent_smiles":"O=C1CCC2OC1OC2","logS":0.27,"solubility_g_L":750.0,"solute_name":"Azoxystrobin","solvent_name":"Cireno (Cyrene)","temperature_K":323.15},
    {"solute_smiles":"OC(C1(CC1)Cl)(CC2=C(Cl)C=CC=C2)CN3N=CNC3=S","solvent_smiles":"OCC1OC(OC1)(C)CC","logS":0.24,"solubility_g_L":600.0,"solute_name":"Prothioconazole","solvent_name":"GBK","temperature_K":298.15},
    {"solute_smiles":"C/C(=N\OCC1=CC=CC=C1/C(=N/OC)/C(=O)OC)/C2=CC(=CC=C2)C(F)(F)F","solvent_smiles":"O=C(OCC)C(=O)C","logS":0.234012307,"solubility_g_L":700.0,"solute_name":"Trifloxystrobin","solvent_name":"Ethyl pyruvate","temperature_K":323.15},
    {"solute_smiles":"OC(C(C)(C)C)(CCC1=CC=C(Cl)C=C1)CN2C=NC=N2","solvent_smiles":"O=C(OCC)CC(=O)C","logS":0.210701389,"solubility_g_L":500.0,"solute_name":"Tebuconazole","solvent_name":"ethyl acetoacetate","temperature_K":323.15},
    {"solute_smiles":"O=C(OC)/C(C1=CC=CC=C1OC2=NC=NC(OC3=CC=CC=C3C#N)=C2)=C/OC","solvent_smiles":"CN1C(CCC1)=O","logS":0.17,"solubility_g_L":590.0,"solute_name":"Azoxystrobin","solvent_name":"NMP","temperature_K":298.15},
    {"solute_smiles":"OC(C(C)(C)C)(CCC1=CC=C(Cl)C=C1)CN2C=NC=N2","solvent_smiles":"O=C(OCC)C(=O)C","logS":0.164943898,"solubility_g_L":450.0,"solute_name":"Tebuconazole","solvent_name":"Ethyl pyruvate","temperature_K":323.15},
    {"solute_smiles":"CNC(=O)C(=NOC)C1=CC=CC=C1OC2=CC=CC=C2","solvent_smiles":"N#CCC(=O)OCC","logS":0.148283132,"solubility_g_L":400.0,"solute_name":"Metominostrobin","solvent_name":"Ethyl cyanoacetate","temperature_K":323.15},
    {"solute_smiles":"O=C1N(N)C(SC)=NN=C1C(C)(C)C","solvent_smiles":"O=C(OCC)C(=O)C","logS":0.15,"solubility_g_L":300.0,"solute_name":"Metribuzin","solvent_name":"Ethyl pyruvate","temperature_K":298.15},
    {"solute_smiles":"CCNC1=NC(NC(C)C)=NC(Cl)=N1","solvent_smiles":"CN1C(CCC1)=O","logS":0.14327111,"solubility_g_L":300.0,"solute_name":"Atrazine","solvent_name":"NMP","temperature_K":323.15},
    {"solute_smiles":"CO/N=C(\C1=CC=CC=C1OC2=C(C(=NC=N2)OC3=CC=CC=C3Cl)F)/C4=NOCCO4","solvent_smiles":"O=C(OCC)C(=O)C","logS":0.12,"solubility_g_L":600.0,"solute_name":"Fluoxastrobin","solvent_name":"Ethyl pyruvate","temperature_K":323.15},
    {"solute_smiles":"C1=CC=C(C(=C1)[C@@H]2[C@@](O2)(CN3C=NC=N3)C4=CC=C(C=C4)F)Cl","solvent_smiles":"O=C(OCC)C(=O)C","logS":0.08,"solubility_g_L":400.0,"solute_name":"Epoxiconazole","solvent_name":"Ethyl pyruvate","temperature_K":323.15},
    {"solute_smiles":"CO/N=C(\C1=CC=CC=C1OC2=C(C(=NC=N2)OC3=CC=CC=C3Cl)F)/C4=NOCCO4","solvent_smiles":"N#CCC(=O)OCC","logS":0.04,"solubility_g_L":500.0,"solute_name":"Fluoxastrobin","solvent_name":"Ethyl cyanoacetate","temperature_K":323.15},
    {"solute_smiles":"CC(C1CC1)C(CN2C=NC=N2)(C3=CC=C(C=C3)Cl)O","solvent_smiles":"N#CCC(=O)OCC","logS":0.01,"solubility_g_L":300.0,"solute_name":"Cyproconazole","solvent_name":"Ethyl cyanoacetate","temperature_K":323.15},
    {"solute_smiles":"CC(C1CC1)C(CN2C=NC=N2)(C3=CC=C(C=C3)Cl)O","solvent_smiles":"O=C(OCC)C(=O)C","logS":0.01,"solubility_g_L":300.0,"solute_name":"Cyproconazole","solvent_name":"Ethyl pyruvate","temperature_K":298.15},
    {"solute_smiles":"OC(C1(CC1)Cl)(CC2=C(Cl)C=CC=C2)CN3N=CNC3=S","solvent_smiles":"O=C(OCC)C(=O)C","logS":0.01,"solubility_g_L":350.0,"solute_name":"Prothioconazole","solvent_name":"Ethyl pyruvate","temperature_K":323.15},
    {"solute_smiles":"C/C(=N\OCC1=CC=CC=C1/C(=N/OC)/C(=O)OC)/C2=CC(=CC=C2)C(F)(F)F","solvent_smiles":"N#CCC(=O)OCC","logS":-0.009025742,"solubility_g_L":400.0,"solute_name":"Trifloxystrobin","solvent_name":"Ethyl cyanoacetate","temperature_K":323.15},
    {"solute_smiles":"O=[N+](/N=C1NCCN\1CC2=CC=C(Cl)N=C2)[O-]","solvent_smiles":"O=C1CCC2OC1OC2","logS":-0.01,"solubility_g_L":250.0,"solute_name":"Imidacloprid","solvent_name":"Cireno (Cyrene)","temperature_K":323.15},
    {"solute_smiles":"OC(C(C)(C)C)(CCC1=CC=C(Cl)C=C1)CN2C=NC=N2","solvent_smiles":"O=C1CCC2OC1OC2","logS":-0.011147361,"solubility_g_L":300.0,"solute_name":"Tebuconazole","solvent_name":"Cireno (Cyrene)","temperature_K":323.15},
    {"solute_smiles":"OC(C(C)(C)C)(CCC1=CC=C(Cl)C=C1)CN2C=NC=N2","solvent_smiles":"N#CCC(=O)OCC","logS":-0.011147361,"solubility_g_L":300.0,"solute_name":"Tebuconazole","solvent_name":"Ethyl cyanoacetate","temperature_K":323.15},
    {"solute_smiles":"OC(C(C)(C)C)(CCC1=CC=C(Cl)C=C1)CN2C=NC=N2","solvent_smiles":"CC(=O)OCC1=CC=CC=C1","logS":-0.011147361,"solubility_g_L":300.0,"solute_name":"Tebuconazole","solvent_name":"Surfonic AG 1705","temperature_K":323.15},
    {"solute_smiles":"OC(C(C)(C)C)(CCC1=CC=C(Cl)C=C1)CN2C=NC=N2","solvent_smiles":"O=C(OCCOC=1C=CC=CC1)C","logS":-0.011147361,"solubility_g_L":300.0,"solute_name":"Tebuconazole","solvent_name":"TEX MEX","temperature_K":323.15},
    {"solute_smiles":"O=C1C(C(C2=CC=C(S(C)(=O)=O)C=C2[N+]([O-])=O)=O)C(CCC1)=O","solvent_smiles":"O=C1CCC2OC1OC2","logS":-0.05,"solubility_g_L":300.0,"solute_name":"Mesotrione","solvent_name":"Cireno (Cyrene)","temperature_K":323.15},
    {"solute_smiles":"OC(C1(CC1)Cl)(CC2=C(Cl)C=CC=C2)CN3N=CNC3=S","solvent_smiles":"O=C1CCC2OC1OC2","logS":-0.06,"solubility_g_L":300.0,"solute_name":"Prothioconazole","solvent_name":"Cireno (Cyrene)","temperature_K":323.15},
    {"solute_smiles":"O=C1C(C(C2=CC=C(S(C)(=O)=O)C=C2[N+]([O-])=O)=O)C(CCC1)=O","solvent_smiles":"CN1C(CCC1)=O","logS":-0.06,"solubility_g_L":295.0,"solute_name":"Mesotrione","solvent_name":"NMP","temperature_K":298.15},
    {"solute_smiles":"O=C1N(CC2=CC=CC=C2Cl)OCC1(C)C","solvent_smiles":"UVCB (Unknown or Variable Composition, Complex Reaction Products)","logS":-0.09,"solubility_g_L":193.33,"solute_name":"Clomazone","solvent_name":"Ultrafilm 2770","temperature_K":298.15},
    {"solute_smiles":"O=C(OC)/C(C1=CC=CC=C1OC2=NC=NC(OC3=CC=CC=C3C#N)=C2)=C/OC","solvent_smiles":"N#CCC(=O)OCC","logS":-0.13,"solubility_g_L":300.0,"solute_name":"Azoxystrobin","solvent_name":"Ethyl cyanoacetate","temperature_K":323.15},
    {"solute_smiles":"O=C(OC)/C(C1=CC=CC=C1OC2=NC=NC(OC3=CC=CC=C3C#N)=C2)=C/OC","solvent_smiles":"O=C(OCC)C(=O)C","logS":-0.13,"solubility_g_L":300.0,"solute_name":"Azoxystrobin","solvent_name":"Ethyl pyruvate","temperature_K":323.15},
    {"solute_smiles":"OC(C1(CC1)Cl)(CC2=C(Cl)C=CC=C2)CN3N=CNC3=S","solvent_smiles":"O=C(OCC)CC(=O)C","logS":-0.14,"solubility_g_L":250.0,"solute_name":"Prothioconazole","solvent_name":"ethyl acetoacetate","temperature_K":323.15},
    {"solute_smiles":"CNC(=O)C(=NOC)C1=CC=CC=C1OC2=CC=CC=C2","solvent_smiles":"N#CCC(=O)OCC","logS":-0.15,"solubility_g_L":200.0,"solute_name":"Metominostrobin","solvent_name":"Ethyl cyanoacetate","temperature_K":298.15},
    {"solute_smiles":"CNC(=O)C(=NOC)C1=CC=CC=C1OC2=CC=CC=C2","solvent_smiles":"O=C(OCC)C(=O)C","logS":-0.152746864,"solubility_g_L":200.0,"solute_name":"Metominostrobin","solvent_name":"Ethyl pyruvate","temperature_K":298.15},
    {"solute_smiles":"OC(C(C)(C)C)(CCC1=CC=C(Cl)C=C1)CN2C=NC=N2","solvent_smiles":"O=C(OCCCC)CCC(=O)C","logS":-0.18723862,"solubility_g_L":200.0,"solute_name":"Tebuconazole","solvent_name":"Butyl Levulinate","temperature_K":298.15},
    {"solute_smiles":"O=C(C1=CC(Br)=NN1C2=NC=CC=C2Cl)NC3=C(C(NC)=O)C=C(Cl)C=C3C","solvent_smiles":"CN1C(CCC1)=O","logS":-0.21,"solubility_g_L":300.0,"solute_name":"Chlorantraniliprole","solvent_name":"NMP","temperature_K":298.15},
    {"solute_smiles":"O=C(C1=CC(Br)=NN1C2=NC=CC=C2Cl)NC3=C(C(NC)=O)C=C(Cl)C=C3C","solvent_smiles":"CN1C(CCC1)=O","logS":-0.21,"solubility_g_L":300.0,"solute_name":"Chlorantraniliprole","solvent_name":"NMP","temperature_K":323.15},
    {"solute_smiles":"O=C(OC)/C(C1=CC=CC=C1OC2=NC=NC(OC3=CC=CC=C3C#N)=C2)=C/OC","solvent_smiles":"O=C(OCC)CC(=O)C","logS":-0.207795885,"solubility_g_L":250.0,"solute_name":"Azoxystrobin","solvent_name":"ethyl acetoacetate","temperature_K":323.15},
    {"solute_smiles":"C1=CC=C(C(=C1)[C@@H]2[C@@](O2)(CN3C=NC=N3)C4=CC=C(C=C4)F)Cl","solvent_smiles":"N#CCC(=O)OCC","logS":-0.22,"solubility_g_L":200.0,"solute_name":"Epoxiconazole","solvent_name":"Ethyl cyanoacetate","temperature_K":323.15},
    {"solute_smiles":"OC(C1(CC1)Cl)(CC2=C(Cl)C=CC=C2)CN3N=CNC3=S","solvent_smiles":"O=C1CCC2OC1OC2","logS":-0.23,"solubility_g_L":201.67,"solute_name":"Prothioconazole","solvent_name":"Cireno (Cyrene)","temperature_K":298.15},
    {"solute_smiles":"OC(C1(CC1)Cl)(CC2=C(Cl)C=CC=C2)CN3N=CNC3=S","solvent_smiles":"CCCCCCc1ccccc1","logS":-0.24,"solubility_g_L":200.0,"solute_name":"Prothioconazole","solvent_name":"Solvesso 150","temperature_K":323.15},
    {"solute_smiles":"OC(C1(CC1)Cl)(CC2=C(Cl)C=CC=C2)CN3N=CNC3=S","solvent_smiles":"CN1C(CCC1)=O","logS":-0.24,"solubility_g_L":198.33,"solute_name":"Prothioconazole","solvent_name":"NMP","temperature_K":298.15},
    {"solute_smiles":"OC(C1(CC1)Cl)(CC2=C(Cl)C=CC=C2)CN3N=CNC3=S","solvent_smiles":"O=C(OCCCC)CCC(=O)C","logS":-0.24,"solubility_g_L":196.67,"solute_name":"Prothioconazole","solvent_name":"Butyl Levulinate","temperature_K":298.15},
    {"solute_smiles":"O=C(OC)/C(C1=CC=CC=C1OC2=NC=NC(OC3=CC=CC=C3C#N)=C2)=C/OC","solvent_smiles":"O=C1CCC2OC1OC2","logS":-0.3,"solubility_g_L":200.0,"solute_name":"Azoxystrobin","solvent_name":"Cireno (Cyrene)","temperature_K":298.15},
    {"solute_smiles":"O=C(OC)/C(C1=CC=CC=C1OC2=NC=NC(OC3=CC=CC=C3C#N)=C2)=C/OC","solvent_smiles":"CC(=O)OCC1=CC=CC=C1","logS":-0.3,"solubility_g_L":200.0,"solute_name":"Azoxystrobin","solvent_name":"Surfonic AG 1705","temperature_K":323.15},
    {"solute_smiles":"OC(C1(CC1)Cl)(CC2=C(Cl)C=CC=C2)CN3N=CNC3=S","solvent_smiles":"CC(=O)OCC1=CC=CC=C1","logS":-0.316175369,"solubility_g_L":166.25,"solute_name":"Prothioconazole","solvent_name":"Surfonic AG 1705","temperature_K":298.15},
    {"solute_smiles":"O=[N+](/N=C1NCCN\1CC2=CC=C(Cl)N=C2)[O-]","solvent_smiles":"O=C1CCC2OC1OC2","logS":-0.33,"solubility_g_L":120.0,"solute_name":"Imidacloprid","solvent_name":"Cireno (Cyrene)","temperature_K":298.15},
    {"solute_smiles":"OC(C1(CC1)Cl)(CC2=C(Cl)C=CC=C2)CN3N=CNC3=S","solvent_smiles":"UVCB (Unknown or Variable Composition, Complex Reaction Products)","logS":-0.33,"solubility_g_L":161.25,"solute_name":"Prothioconazole","solvent_name":"Ultrafilm 2770","temperature_K":298.15},
    {"solute_smiles":"O=C1N(N)C(SC)=NN=C1C(C)(C)C","solvent_smiles":"N#CCC(=O)OCC","logS":-0.33,"solubility_g_L":100.0,"solute_name":"Metribuzin","solvent_name":"Ethyl cyanoacetate","temperature_K":298.15},
    {"solute_smiles":"CCNC1=NC(NC(C)C)=NC(Cl)=N1","solvent_smiles":"O=C(OCC)C(=O)C","logS":-0.333850145,"solubility_g_L":100.0,"solute_name":"Atrazine","solvent_name":"Ethyl pyruvate","temperature_K":323.15},
    {"solute_smiles":"O=C1C(C(C2=CC=C(S(C)(=O)=O)C=C2[N+]([O-])=O)=O)C(CCC1)=O","solvent_smiles":"O=C1CCC2OC1OC2","logS":-0.35,"solubility_g_L":153.0,"solute_name":"Mesotrione","solvent_name":"Cireno (Cyrene)","temperature_K":298.15},
    {"solute_smiles":"OC(C1(CC1)Cl)(CC2=C(Cl)C=CC=C2)CN3N=CNC3=S","solvent_smiles":"O=C(OCC)CC(=O)C","logS":-0.36,"solubility_g_L":150.0,"solute_name":"Prothioconazole","solvent_name":"ethyl acetoacetate","temperature_K":298.15},
    {"solute_smiles":"OC(C1(CC1)Cl)(CC2=C(Cl)C=CC=C2)CN3N=CNC3=S","solvent_smiles":"N#CCC(=O)OCC","logS":-0.36,"solubility_g_L":150.0,"solute_name":"Prothioconazole","solvent_name":"Ethyl cyanoacetate","temperature_K":323.15},
    {"solute_smiles":"CCC(C)C1C(C=CC2(O1)CC3CC(O2)CC=C(C(C(C=CC=C4COC5C4(C(C=C(C5O)C)C(=O)O3)O)C)OC6CC(C(C(O6)C)OC7CC(C(C(O7)C)O)OC)OC)C)C","solvent_smiles":"O=C(OCC)C(=O)C","logS":-0.396995944,"solubility_g_L":350.0,"solute_name":"Abamectin","solvent_name":"Ethyl pyruvate","temperature_K":323.15},
    {"solute_smiles":"O=[N+](/N=C1NCCN\1CC2=CC=C(Cl)N=C2)[O-]","solvent_smiles":"N#CCC(=O)OCC","logS":-0.41,"solubility_g_L":100.0,"solute_name":"Imidacloprid","solvent_name":"Ethyl cyanoacetate","temperature_K":323.15},
    {"solute_smiles":"OC(C1(CC1)Cl)(CC2=C(Cl)C=CC=C2)CN3N=CNC3=S","solvent_smiles":"O=C(OCC)C=1C=CC=CC1","logS":-0.47,"solubility_g_L":118.0,"solute_name":"Prothioconazole","solvent_name":"Ethyl Benzoate","temperature_K":298.15},
    {"solute_smiles":"CC(C1CC1)C(CN2C=NC=N2)(C3=CC=C(C=C3)Cl)O","solvent_smiles":"N#CCC(=O)OCC","logS":-0.47,"solubility_g_L":100.0,"solute_name":"Cyproconazole","solvent_name":"Ethyl cyanoacetate","temperature_K":298.15},
    {"solute_smiles":"OC(C(C)(C)C)(CCC1=CC=C(Cl)C=C1)CN2C=NC=N2","solvent_smiles":"O=C(OCC)CC(=O)C","logS":-0.488268615,"solubility_g_L":100.0,"solute_name":"Tebuconazole","solvent_name":"ethyl acetoacetate","temperature_K":298.15},
    {"solute_smiles":"C1=CC=C(C(=C1)[C@@H]2[C@@](O2)(CN3C=NC=N3)C4=CC=C(C=C4)F)Cl","solvent_smiles":"O=C(OCC)C(=O)C","logS":-0.52,"solubility_g_L":100.0,"solute_name":"Epoxiconazole","solvent_name":"Ethyl pyruvate","temperature_K":298.15},
    {"solute_smiles":"O=C1C(C(C2=CC=C(S(C)(=O)=O)C=C2[N+]([O-])=O)=O)C(CCC1)=O","solvent_smiles":"N#CCC(=O)OCC","logS":-0.53,"solubility_g_L":100.0,"solute_name":"Mesotrione","solvent_name":"Ethyl cyanoacetate","temperature_K":323.15},
    {"solute_smiles":"O=C1C(C(C2=CC=C(S(C)(=O)=O)C=C2[N+]([O-])=O)=O)C(CCC1)=O","solvent_smiles":"O=C(OCC)C(=O)C","logS":-0.53,"solubility_g_L":100.0,"solute_name":"Mesotrione","solvent_name":"Ethyl pyruvate","temperature_K":323.15},
    {"solute_smiles":"OC(C1(CC1)Cl)(CC2=C(Cl)C=CC=C2)CN3N=CNC3=S","solvent_smiles":"N#CCC(=O)OCC","logS":-0.54,"solubility_g_L":100.0,"solute_name":"Prothioconazole","solvent_name":"Ethyl cyanoacetate","temperature_K":298.15},
    {"solute_smiles":"OC(C1(CC1)Cl)(CC2=C(Cl)C=CC=C2)CN3N=CNC3=S","solvent_smiles":"O=C(OCC)C(=O)C","logS":-0.54,"solubility_g_L":100.0,"solute_name":"Prothioconazole","solvent_name":"Ethyl pyruvate","temperature_K":298.15},
    {"solute_smiles":"OC(C1(CC1)Cl)(CC2=C(Cl)C=CC=C2)CN3N=CNC3=S","solvent_smiles":"O=C(OCCOC=1C=CC=CC1)C","logS":-0.54,"solubility_g_L":100.0,"solute_name":"Prothioconazole","solvent_name":"TEX MEX","temperature_K":323.15},
    {"solute_smiles":"O=C(OC)/C(C1=CC=CC=C1OC2=NC=NC(OC3=CC=CC=C3C#N)=C2)=C/OC","solvent_smiles":"O=C(OCCOC=1C=CC=CC1)C","logS":-0.61,"solubility_g_L":100.0,"solute_name":"Azoxystrobin","solvent_name":"TEX MEX","temperature_K":323.15},
    {"solute_smiles":"C/C(=N\OCC1=CC=CC=C1/C(=N/OC)/C(=O)OC)/C2=CC(=CC=C2)C(F)(F)F","solvent_smiles":"N#CCC(=O)OCC","logS":-0.611085733,"solubility_g_L":100.0,"solute_name":"Trifloxystrobin","solvent_name":"Ethyl cyanoacetate","temperature_K":298.15},
    {"solute_smiles":"C/C(=N\OCC1=CC=CC=C1/C(=N/OC)/C(=O)OC)/C2=CC(=CC=C2)C(F)(F)F","solvent_smiles":"O=C(OCC)C(=O)C","logS":-0.611085733,"solubility_g_L":100.0,"solute_name":"Trifloxystrobin","solvent_name":"Ethyl pyruvate","temperature_K":298.15},
    {"solute_smiles":"CCC(C)C1C(C=CC2(O1)CC3CC(O2)CC=C(C(C(C=CC=C4COC5C4(C(C=C(C5O)C)C(=O)O3)O)C)OC6CC(C(C(O6)C)OC7CC(C(C(O7)C)O)OC)OC)C)C","solvent_smiles":"CN1C(CCC1)=O","logS":-0.63292661,"solubility_g_L":203.3,"solute_name":"Abamectin","solvent_name":"NMP","temperature_K":298.15},
    {"solute_smiles":"CCC(C)C1C(C=CC2(O1)CC3CC(O2)CC=C(C(C(C=CC=C4COC5C4(C(C=C(C5O)C)C(=O)O3)O)C)OC6CC(C(C(O6)C)OC7CC(C(C(O7)C)O)OC)OC)C)C","solvent_smiles":"CC(=O)OCC1=CC=CC=C1","logS":-0.640033993,"solubility_g_L":200.0,"solute_name":"Abamectin","solvent_name":"Surfonic AG 1705","temperature_K":323.15},
    {"solute_smiles":"CO/N=C(\C1=CC=CC=C1OC2=C(C(=NC=N2)OC3=CC=CC=C3Cl)F)/C4=NOCCO4","solvent_smiles":"N#CCC(=O)OCC","logS":-0.66,"solubility_g_L":100.0,"solute_name":"Fluoxastrobin","solvent_name":"Ethyl cyanoacetate","temperature_K":298.15},
    {"solute_smiles":"CO/N=C(\C1=CC=CC=C1OC2=C(C(=NC=N2)OC3=CC=CC=C3Cl)F)/C4=NOCCO4","solvent_smiles":"O=C(OCC)C(=O)C","logS":-0.66,"solubility_g_L":100.0,"solute_name":"Fluoxastrobin","solvent_name":"Ethyl pyruvate","temperature_K":298.15},
    {"solute_smiles":"FC1=CC=CC(F)=C1C(NC(NC2=C(Cl)C=C(OC(F)(F)C(C(F)(F)F)F)C(Cl)=C2)=O)=O","solvent_smiles":"O=C(OCC)C(=O)C","logS":-0.71,"solubility_g_L":100.0,"solute_name":"Lufenuron","solvent_name":"Ethyl pyruvate","temperature_K":323.15},
    {"solute_smiles":"CCC(C)C1C(C=CC2(O1)CC3CC(O2)CC=C(C(C(C=CC=C4COC5C4(C(C=C(C5O)C)C(=O)O3)O)C)OC6CC(C(C(O6)C)OC7CC(C(C(O7)C)O)OC)OC)C)C","solvent_smiles":"OCC1OC(OC1)(C)CC","logS":-0.757794145,"solubility_g_L":152.5,"solute_name":"Abamectin","solvent_name":"GBK","temperature_K":298.15},
    {"solute_smiles":"CCC(C)C1C(C=CC2(O1)CC3CC(O2)CC=C(C(C(C=CC=C4COC5C4(C(C=C(C5O)C)C(=O)O3)O)C)OC6CC(C(C(O6)C)OC7CC(C(C(O7)C)O)OC)OC)C)C","solvent_smiles":"CC(=O)OCC1=CC=CC=C1","logS":-0.872878126,"solubility_g_L":117.0,"solute_name":"Abamectin","solvent_name":"Surfonic AG 1705","temperature_K":298.15},
    {"solute_smiles":"CCC(C)C1C(C=CC2(O1)CC3CC(O2)CC=C(C(C(C=CC=C4COC5C4(C(C=C(C5O)C)C(=O)O3)O)C)OC6CC(C(C(O6)C)OC7CC(C(C(O7)C)O)OC)OC)C)C","solvent_smiles":"O=C1CCC2OC1OC2","logS":-0.941063988,"solubility_g_L":100.0,"solute_name":"Abamectin","solvent_name":"Cireno (Cyrene)","temperature_K":323.15},
])

])

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

