import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import numpy as np

# === Load Input Data (Equivalent to MATLAB readtable) === #
path_main = "C:\\Users\\nanf\\Desktop\\H2-districts\\Input"


# Load Excel & CSV files
Input1 = pd.read_excel(f"{path_main}/Inputs_BoundaryLoads_2022.xlsx")
Input2 = pd.read_csv(f"{path_main}/nestel_demand.csv")
Input3 = pd.read_excel(f"{path_main}/PriceImpExp_ele.xlsx")

# Extract relevant data
irradiance = Input1["G"].values
P_load = Input2["P_el_nest_h"].values
c_gridimp = Input3["Priceimpgreen"].values
c_gridexp = Input3["PriceexpPV"].values

nHours = 8760  # Total hours simulated

# === Problem Parameters === #
bigM = 600
MaxSimTime = 1000  # Max time for solver (seconds)
deltat = 3600  # Time step = 1h (in seconds)

HHV = 39.39 * 3.6 * 10 ** 3  # Hydrogen higher heating value (kJ/kg)
c_h2 = 4  # Hydrogen cost (CHF/kg)
eff_b = 0.95  # Battery efficiency
eff_PV = 0.17  # PV efficiency
eff_FC = 0.62  # Fuel Cell efficiency

nH_lifetime = 40000  # PEMFC lifetime (hours)
Delta_eff_FC = 0.1 * eff_FC
Const1 = Delta_eff_FC / nH_lifetime
Const3 = Const1 / (HHV * eff_FC ** 2)

# Sizing Limits
S_FC_min, S_FC_max = 0, 60  # Fuel Cell size range (kW)
P_PV_peak_ref = 100  # Peak PV power (kW)
Area_PV_ref = P_PV_peak_ref / eff_PV
Area_PV_min, Area_PV_max = 0, 2 * Area_PV_ref  # PV area range (m²)
C_b_max = 96 * 3.6 * 10 ** 3  # Battery capacity (kJ)
C_b_min = C_b_max / 3  # Min Battery capacity
P_b_max = C_b_max / 3600  # Max battery power (kW)

# === Cost Parameters === #
d = 0.05
ann = d / (1 - (1 + d) ** -20)  # Annuity factor

UP_PV, life_PV, maint_PV = 1300, 25, 0.01
ann_PV = d / (1 - (1 + d) ** -life_PV)
UP_b, life_b, maint_b = 1000, 12, 0.02
ann_b = d / (1 - (1 + d) ** -life_b)
UP_FC, life_FC, maint_FC = 950, 7.5, 0.024
ann_FC = d / (1 - (1 + d) ** -life_FC)

# === Create Gurobi Model === #
model = gp.Model("MIQP_Energy_System")

# === Decision Variables === #
C_b = model.addVar(lb=C_b_min, ub=C_b_max, vtype=GRB.CONTINUOUS, name="C_b")  # Battery capacity (kJ)
S_FC = model.addVar(lb=S_FC_min, ub=S_FC_max, vtype=GRB.CONTINUOUS, name="S_FC")  # Fuel Cell size (kW)
Area_PV = model.addVar(lb=Area_PV_min, ub=Area_PV_max, vtype=GRB.CONTINUOUS, name="Area_PV")  # PV area (m²)

# Operational Variables
P_FC = model.addVars(nHours, lb=0, ub=S_FC_max, vtype=GRB.CONTINUOUS, name="P_FC")
P_FC_On = model.addVars(nHours, lb=0, ub=S_FC_max, vtype=GRB.CONTINUOUS, name="P_FC_On")
FC_On = model.addVars(nHours, vtype=GRB.BINARY, name="FC_On")

w_hours = model.addVars(nHours, lb=0, ub=8760, vtype=GRB.CONTINUOUS, name="w_hours")

P_imp = model.addVars(nHours, lb=0, vtype=GRB.CONTINUOUS, name="P_imp")
P_exp = model.addVars(nHours, lb=0, vtype=GRB.CONTINUOUS, name="P_exp")

E_b = model.addVars(nHours, lb=0, ub=C_b_max, vtype=GRB.CONTINUOUS, name="E_b")
P_b_ch = model.addVars(nHours, lb=0, ub=P_b_max, vtype=GRB.CONTINUOUS, name="P_b_ch")
P_b_disch = model.addVars(nHours, lb=0, ub=P_b_max, vtype=GRB.CONTINUOUS, name="P_b_disch")

# Derived Variables
P_PV = {t: irradiance[t] * eff_PV * Area_PV / 1000 for t in range(nHours)}

# === Constraints === #
for t in range(nHours):
    model.addConstr(P_FC[t] + P_PV[t] + P_b_disch[t] + P_imp[t] == P_b_ch[t] + P_load[t] + P_exp[t])

    if t > 0:
        model.addConstr(E_b[t] == E_b[t - 1] + P_b_ch[t] * deltat - P_b_disch[t] * deltat)

    model.addConstr(P_FC[t] <= P_FC_On[t])
    model.addConstr(P_FC_On[t] <= S_FC_max * FC_On[t])
    model.addConstr(P_FC_On[t] <= S_FC)

    if t > 0:
        model.addConstr(w_hours[t] == w_hours[t - 1] + FC_On[t])

# === Quadratic Objective Function === #
m_flow_H2 = {t: P_FC[t] / eff_FC / HHV * 3600 + P_FC[t] * Const3 * w_hours[t] * 3600 for t in range(nHours)}

cost_inst = (P_PV_peak_ref * UP_PV * ann_PV + S_FC * UP_FC * ann_FC + C_b / 3600 * UP_b * ann_b) / 1000
cost_imp = gp.quicksum(c_h2 * m_flow_H2[t] for t in range(nHours)) / 1000 + gp.quicksum(
    P_imp[t] * c_gridimp[t] for t in range(nHours)) / 1000
cost_exp = gp.quicksum(P_exp[t] * c_gridexp[t] for t in range(nHours)) / 1000
cost_maint = (maint_PV * P_PV_peak_ref * UP_PV + maint_FC * S_FC * UP_FC + maint_b * UP_b * C_b / 3600) / 1000

model.setObjective(cost_inst + cost_imp - cost_exp + cost_maint, GRB.MINIMIZE)

# === Solve the MIQP Problem === #
model.optimize()

# === Display Results === #
print("Optimal PV Area:", Area_PV.X)
print("Optimal Fuel Cell Size:", S_FC.X)
print("Optimal Battery Capacity:", C_b.X)
