"""Track refitted equilibria as power-carrying flow is introduced.

Alpha is a numerical continuation coordinate, held fixed at each refit. It is
not an additional calibration parameter. The final closure uses alpha=1.
"""

import argparse
from dataclasses import asdict, replace
import json
import time

import numpy as np
from scipy.optimize import least_squares

from flow_consistency import Experiment, HERE, fit


def run(baseline=False, corrected_drive=True, stage2_seed=False, loss_seed=False):
    ex = Experiment(baseline, corrected_drive)
    plates = np.array([r["stage_temperature_K"] for r in ex.saved["points"]])
    loads = np.array([r["heater_power_W"] for r in ex.saved["points"]])
    states = np.array([ex.state(i) for i in range(6)])
    p0 = ex.parameters
    if stage2_seed:
        saved = json.loads((HERE/"fit_2.json").read_text())
        p0 = fit.Parameters(**saved["parameters"])
        states = np.array([np.r_[np.log(np.array(r["temperature_K"])[ex.nodes]),
            r["mean_pressure_Pa"]/1e6,r["swing_Pa"]/1e6] for r in saved["points"]])
    if loss_seed:
        saved = json.loads((HERE/"quadrature_baseline_1.json").read_text())
        p0 = fit.Parameters(**saved["parameters"])
        states = np.array([r["state"] for r in saved["points"]])
    n = states.shape[1]
    coefficients = np.r_[p0.packing_constant,np.array(p0.conductance)/1e-10]
    if not baseline:
        coefficients = np.r_[coefficients,p0.exchanger_effectiveness]

    def parameters(c):
        return replace(p0,packing_constant=tuple(c[:2]),conductance=tuple(c[2:4]*1e-10),
            exchanger_effectiveness=p0.exchanger_effectiveness if baseline else c[4])

    cap = min(ex.drive.relief_swing)/1e6
    lower_state = np.r_[np.full(n-2,np.log(2)),cap/2+.01,.001]
    upper_state = np.r_[np.full(n-2,np.log(300)),2.5,cap]
    lower = np.r_[np.tile(lower_state,6),[1e-7,1e-7,0,1e-7],[] if baseline else [1e-5]]
    upper = np.r_[np.tile(upper_state,6),[10,10,100,100],[] if baseline else [100]]
    x = np.clip(np.r_[states.ravel(),coefficients],lower+1e-10,upper-1e-10)
    started = time.perf_counter()
    history = []
    last_alpha = 0
    for alpha in ((1,) if loss_seed else (0,.001,.01,.03,.1,.3,1)):
        def point_residual(state,c,index):
            balance,values,_ = ex.evaluate(state,plates[index],parameters(c),alpha)
            return np.r_[balance*1e4,values[:2]-loads[index]]

        def residual(candidate):
            return np.concatenate([point_residual(s,candidate[6*n:],i)
                for i,s in enumerate(candidate[:6*n].reshape(6,n))])

        def jacobian(candidate):
            jac = np.zeros((6*(n+2),len(candidate)))
            c = candidate[6*n:]
            for i,s in enumerate(candidate[:6*n].reshape(6,n)):
                rows = slice(i*(n+2),(i+1)*(n+2))
                for k in range(n):
                    delta = np.zeros(n)
                    delta[k] = 2e-6*max(abs(s[k]),1)
                    jac[rows,i*n+k] = (point_residual(s+delta,c,i)-point_residual(s-delta,c,i))/(2*delta[k])
                for k in range(len(c)):
                    delta = np.zeros(len(c))
                    delta[k] = 2e-6*max(abs(c[k]),.01)
                    jac[rows,6*n+k] = (point_residual(s,c+delta,i)-point_residual(s,c-delta,i))/(2*delta[k])
            return jac

        result = least_squares(residual,x,jac=jacobian,bounds=(lower,upper),x_scale="jac",
            max_nfev=80,ftol=1e-9,xtol=1e-9,gtol=1e-8)
        p = parameters(result.x[6*n:])
        candidate_states,errors,balances = [],[],[]
        failure = None
        for i,s in enumerate(result.x[:6*n].reshape(6,n)):
            try:
                state,values,details = ex.solve(s,plates[i],p,alpha)
            except RuntimeError as error:
                failure = str(error)
                break
            candidate_states.append(state)
            errors.append(values[:2]-loads[i])
            balances.append(float(max(abs(details["heat"][ex.nodes]))))
        row = dict(alpha=alpha,parameters=asdict(p),optimizer_success=result.success,
            nfev=result.nfev,elapsed_s=time.perf_counter()-started,failure=failure)
        approximate = residual(result.x).reshape(6,n+2)
        row["approximate_max_balance_W"] = float(np.max(abs(approximate[:,:n-2]))/100)
        row["minimum_matrix_K"] = np.exp(result.x[:6*n].reshape(6,n)[:,:-2]).min(axis=1).tolist()
        if failure is None:
            row.update(rms_W=np.sqrt(np.mean(np.array(errors)**2,axis=0)).tolist(),
                       max_balance_W=max(balances),residual_W=np.array(errors).tolist())
            x = np.r_[np.array(candidate_states).ravel(),result.x[6*n:]]
            last_alpha = alpha
        history.append(row)
        print(json.dumps(row),flush=True)
        if failure is not None:
            row["failed_states"] = result.x[:6*n].reshape(6,n).tolist()
            break

    p = parameters(x[6*n:])
    points = []
    for i,s in enumerate(x[:6*n].reshape(6,n)):
        balance,values,details = ex.evaluate(s,plates[i],p,last_alpha)
        # Eliminate the two instantaneous drive variables from the thermal Jacobian.
        derivatives = []
        for direction in np.eye(n):
            plus = ex.evaluate(s+1e-5*direction,plates[i],p,last_alpha)[0]
            minus = ex.evaluate(s-1e-5*direction,plates[i],p,last_alpha)[0]
            derivatives.append((plus-minus)/2e-5)
        jac = np.array(derivatives).T
        thermal = jac[:-2,:-2]-jac[:-2,-2:]@np.linalg.solve(jac[-2:,-2:],jac[-2:,:-2])
        capacities = np.concatenate(ex.head.calc_heat_capacity(details["temperature"]))
        growth = np.linalg.eigvals(100*thermal/(capacities*details["temperature"][ex.nodes])[:,None]).real.max()
        point = dict(point=i,state=s.tolist(),temperature_K=details["temperature"].tolist(),
            equilibrium_verified=bool(max(abs(balance))<1e-7),
            residual_W=(values[:2]-loads[i]).tolist(),growth_rate_per_s=float(growth),
            diagnostics=values.tolist())
        if i in (1,4):
            point["sensitivity"] = ex.sensitivity(s,plates[i],alpha=last_alpha,p=p)
        points.append(point)
    output = dict(baseline=baseline,corrected_drive=corrected_drive,stage2_seed=stage2_seed,
        final_alpha=last_alpha,parameters=asdict(p),history=history,points=points)
    output["equilibrium_verified"] = all(r["equilibrium_verified"] for r in points)
    name = f"continuation_{'baseline' if baseline else 'current'}_{'drive' if corrected_drive else 'loss'}"
    if stage2_seed:
        name += "_stage2_seed"
    if loss_seed:
        name += "_loss_seed"
    (HERE/f"{name}.json").write_text(json.dumps(output,indent=2)+"\n")
    print("FINAL",name,"alpha",last_alpha,"slopes",
        [r["sensitivity"]["estimates"][-1]["dT1_dT2_at_fixed_Q1"] for r in points if "sensitivity" in r],flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline",action="store_true")
    parser.add_argument("--loss-only",action="store_true")
    parser.add_argument("--seed-stage2",action="store_true")
    parser.add_argument("--seed-loss-fit",action="store_true")
    args = parser.parse_args()
    run(args.baseline,not args.loss_only,args.seed_stage2,args.seed_loss_fit)
