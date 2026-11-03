"""Audit flow amplitudes and test a quadrature closure in an isolated experiment.

The closure retains the existing drive and compliance proxy. In-phase amplitude
is 2*rho*E/Pa; its equivalent positive-transfer-per-cycle rate is amplitude/pi.
Stage 1 uses constant in-phase mass flow, anchored at the branch to both stages.
This is a reduced diagnostic, not a pressure/flow network.
"""

import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import fit
from pulse_tube_cooldown.library.links.cooler.head import calc_acoustic_power


class Experiment:
    def __init__(self, baseline=False, corrected_drive=False):
        inputs, self.data, self.warm, self.drive = fit.read_experiment()
        self.baseline = baseline
        self.corrected_drive = corrected_drive
        directory = fit.HERE / "baseline" if baseline else fit.HERE
        self.saved = json.loads((directory / "results.json").read_text())
        self.head = fit.read_head(inputs["head"], directory=directory / "head_model" if baseline else None)
        self.nodes = fit.column_indices(self.head)
        self.walls = fit.build_walls()
        self.parameters = fit.Parameters(
            packing_constant=tuple(self.saved["packing_constants"]),
            conductance=tuple(self.saved["conductances_W_per_Pa2"]),
            exchanger_effectiveness=self.saved.get("exchanger_effectiveness", 0.25),
            frequency=self.saved["frequency_Hz"],
        )
        self.drive = replace(self.drive, charge=self.saved["helium_charge_kg"])

    def state(self, point):
        row = self.saved["points"][point]
        t = np.array(row["temperature_K"])
        high, low = row["supply_pressure_Pa"], row["return_pressure_Pa"]
        return np.r_[np.log(t[self.nodes]), (high+low)/2e6, (high-low)/1e6]

    def evaluate(self, state, plate, p=None, alpha=0, mode="full", reference=None):
        p = p or self.parameters
        h, d = self.head, self.drive
        t = np.empty(h.stage_2+1)
        t[[0, h.stage_1, -1]] = [self.warm, *plate]
        t[self.nodes] = np.exp(state[:-2])
        mean, swing = state[-2:]*1e6
        pa = d.waveform_factor*swing/2
        capacity = h.calc_heat_capacity(t)
        eps = fit.calc_capacity_ratio(t, capacity, h.pore_volume, mean)
        gross = tuple(fit.calc_enthalpy_flow(c, swing, d.waveform_factor) for c in p.conductance)
        proxy = fit.calc_mass_flow(t, h.pore_volume, h.pulse_tube_volume, mean+swing/2, mean-swing/2, p.frequency)
        if self.baseline:
            boundary = (t[h.stage_1:-1]+t[h.stage_1+1:])/2
            boundary[[0, -1]] = t[[h.stage_1, -1]]
            rho = fit.helium.density(boundary, mean)
            acoustic = gross[1]*rho[-1]/rho
            endpoint = t[-1]
        else:
            boundary = (t[h.stage_1:-1]+t[h.stage_1+1:])/2
            boundary[[0, -1]] = t[[h.stage_1, -2]]
            rho = fit.helium.density(boundary, mean)
            acoustic = calc_acoustic_power(t, eps, gross[1], mean)
            endpoint = t[-2]
        inphase2 = 2*rho*acoustic/pa
        inphase1 = 2*float(fit.helium.density(t[h.stage_1], mean))*(gross[0]+acoustic[0])/pa
        inphase = np.r_[np.full(h.stage_1, inphase1), inphase2]
        thermal = np.sqrt(proxy**2+alpha*(inphase/np.pi)**2)
        if mode == "normalized_first_flow":
            thermal[:h.stage_1] = reference["proxy"][:h.stage_1]*inphase1/reference["inphase1"]
        if self.baseline:
            zeros = tuple(np.zeros_like(e) for e in eps)
            heat = fit.calc_uncoupled_heat_flow(t, capacity, zeros, thermal, gross,
                p.packing_constant, p.exchanger_effectiveness, mean, p.frequency)
        else:
            heat = fit.calc_heat_flow(t, capacity, eps, thermal, gross,
                p.packing_constant, p.exchanger_effectiveness, mean, p.frequency)
        inventory = fit.calc_gas_mass(t, h.pore_volume, h.pulse_tube_volume, mean)
        inventory += d.warm_volume*float(fit.helium.density(self.warm, mean))
        compressor = float(fit.helium.density(self.warm, mean-swing/2))*sum(d.displacement_rate)
        pressure_residual = (inventory-d.charge)/d.charge
        demand = thermal[0] if self.corrected_drive else proxy[0]
        drive_residual = min((min(d.relief_swing)-swing)/min(d.relief_swing), (compressor-demand)/0.005)
        if mode == "fixed_pressures":
            pressure_residual, drive_residual = (state[-2:]-reference["state"][-2:])
        wall = fit.calc_link_heat_flow(self.walls,
            {fit.WARM_FLANGE:self.warm, fit.STAGE_1:plate[0], fit.STAGE_2:plate[1]})
        cooling = -heat[[h.stage_1, -1]]-[wall[fit.STAGE_1], wall[fit.STAGE_2]]
        k = h.stage_1-1
        cp = float(fit.helium.isobaric_specific_heat((t[k]+t[k+1])/2, mean))
        loss = 2*p.packing_constant[0]*(thermal[k]*cp)**2/(p.frequency*capacity[0][-1])*(t[k]-t[k+1])
        noeps = 2*gross[1]*float(fit.helium.density(endpoint, mean))/pa
        values = np.r_[cooling, sum(gross), loss, inphase1, noeps, endpoint, eps[1][-1],
                       proxy[-1]/p.frequency, inphase2[-1]/(np.pi*p.frequency),
                       compressor/p.frequency, thermal[0]/p.frequency]
        details = dict(proxy=proxy, inphase1=inphase1, state=state, temperature=t,
                       heat=heat, inphase=inphase, gross=gross, acoustic=acoustic)
        return np.r_[heat[self.nodes]/100, pressure_residual, drive_residual], values, details

    def solve(self, state, plate, p=None, alpha=0):
        n = len(self.nodes)
        cap = min(self.drive.relief_swing)/1e6
        lower = np.r_[np.full(n, np.log(2.0)), cap/2+0.01, 0.001]
        upper = np.r_[np.full(n, np.log(300.0)), 2.5, cap]
        solution = least_squares(lambda y:self.evaluate(y,plate,p,alpha)[0],
            np.clip(state,lower+1e-12,upper-1e-12), bounds=(lower,upper),
            x_scale="jac", jac="3-point", diff_step=1e-6,
            xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=200)
        residual, values, details = self.evaluate(solution.x,plate,p,alpha)
        if max(abs(residual)) > 1e-7:
            raise RuntimeError(f"Equilibrium residual {max(abs(residual))}, alpha={alpha}")
        return solution.x, values, details

    def sensitivity(self, state, plate, mode="full", alpha=0, p=None):
        base = self.evaluate(state,plate,p,alpha)
        estimates = []
        for step in (1e-5, 5e-6):
            jac, response = [], []
            for direction in np.eye(len(state)):
                plus = self.evaluate(state+step*direction,plate,p,alpha,mode,base[2])
                minus = self.evaluate(state-step*direction,plate,p,alpha,mode,base[2])
                jac.append((plus[0]-minus[0])/(2*step))
                response.append((plus[1]-minus[1])/(2*step))
            derivatives = []
            for direction in np.eye(2):
                plus = self.evaluate(state,plate+.001*direction,p,alpha,mode,base[2])
                minus = self.evaluate(state,plate-.001*direction,p,alpha,mode,base[2])
                dy = np.linalg.solve(np.array(jac).T,-(plus[0]-minus[0])/.002)
                derivatives.append(np.array(response).T@dy+(plus[1]-minus[1])/.002)
            first, second = derivatives
            estimates.append(dict(step=step,d_dT1=first.tolist(), d_dT2=second.tolist(),
                dT1_dT2_at_fixed_Q1=float(-second[0]/first[0]),
                dQsum_dT2=float(second[0]+second[1])))
        return dict(reference=base[1].tolist(), mode=mode, estimates=estimates)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quadrature", action="store_true")
    args = parser.parse_args()
    output = dict(quantity_order=["Q1_W","Q2_W","Hsum_W","L1_W","inphase1_kg_s",
        "inphase2_noeps_kg_s","endpoint_K","eps_cold","compliance_cold_kg_per_cycle",
        "inphase_cold_kg_per_cycle","compressor_kg_per_cycle","quadrature_warm_kg_per_cycle"], rows=[])
    for baseline in (True,False):
        ex = Experiment(baseline)
        label = "baseline" if baseline else "current"
        for index in (1,4):
            plate = np.array(ex.saved["points"][index]["stage_temperature_K"])
            state = ex.state(index)
            residual, values, details = ex.evaluate(state,plate)
            np.testing.assert_allclose(residual,0,atol=1e-8)
            expected = np.array(ex.saved["points"][index]["heater_power_W"])-ex.saved["points"][index]["residual_W"]
            np.testing.assert_allclose(values[:2],expected,atol=1e-7)
            for mode in ("full","normalized_first_flow","fixed_pressures"):
                row = dict(model=label,point=index,**ex.sensitivity(state,plate,mode))
                output["rows"].append(row)
                print(label,index,mode,row["estimates"][-1],flush=True)
            if args.quadrature:
                continuation = []
                for alpha in np.linspace(.01,1,21):
                    try:
                        state, values, details = ex.solve(state,plate,alpha=alpha)
                    except RuntimeError as error:
                        continuation.append(dict(alpha=float(alpha), error=str(error)))
                        break
                    continuation.append(dict(alpha=float(alpha),values=values.tolist()))
                output["rows"].append(dict(model=label,point=index,quadrature_fixed_parameters=continuation))
                print(label,index,"quadrature",continuation[-1],flush=True)
    output["empirical_group_slopes"] = []
    for group in ([2,1,0],[5,4,3]):
        t = np.array([ex.saved["points"][j]["stage_temperature_K"] for j in group])
        output["empirical_group_slopes"].append(dict(Q1_W=ex.saved["points"][group[0]]["heater_power_W"][0],
            secant_dT1_dT2=float((t[-1,0]-t[0,0])/(t[-1,1]-t[0,1])),
            regression_dT1_dT2=float(np.polyfit(t[:,1],t[:,0],1)[0])))
    (HERE/"flow_consistency.json").write_text(json.dumps(output,indent=2)+"\n")


if __name__ == "__main__":
    main()
