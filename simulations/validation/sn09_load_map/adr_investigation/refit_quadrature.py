"""Simultaneous-state exploratory refit with the quadrature flow diagnostic."""

import argparse
from dataclasses import asdict, replace
import json
import time

import numpy as np
from scipy.optimize import least_squares

from flow_consistency import Experiment, HERE


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline",action="store_true")
    parser.add_argument("--alpha",type=float,default=1)
    parser.add_argument("--max-nfev",type=int,default=150)
    parser.add_argument("--seed-stage2",action="store_true")
    args = parser.parse_args()
    ex = Experiment(args.baseline)
    states = np.array([ex.state(i) for i in range(6)])
    plates = np.array([row["stage_temperature_K"] for row in ex.saved["points"]])
    heater = np.array([row["heater_power_W"] for row in ex.saved["points"]])
    p0 = ex.parameters
    if args.seed_stage2:
        seed = json.loads((HERE/"fit_2.json").read_text())
        p0 = replace(p0,**seed["parameters"])
        states = np.array([np.r_[np.log(np.array(r["temperature_K"])[ex.nodes]),
            r["mean_pressure_Pa"]/1e6,r["swing_Pa"]/1e6] for r in seed["points"]])
    n = states.shape[1]
    gains = []
    for i in range(6):
        details = ex.evaluate(states[i],plates[i],p0)[2]
        gains.append(1+args.alpha*(details["inphase"]/(np.pi*details["proxy"]))**2)
    gains = np.mean(gains,axis=0)
    h = ex.head
    coefficients = np.array([p0.packing_constant[0]/np.mean(gains[:h.stage_1]),
        p0.packing_constant[1]/np.mean(gains[h.stage_1:]),
        p0.conductance[0]/1e-10,p0.conductance[1]/1e-10,
        p0.exchanger_effectiveness/np.sqrt(gains[-1])])
    if args.baseline:
        coefficients = coefficients[:4]

    def parameters(c):
        return replace(p0,packing_constant=tuple(c[:2]),conductance=tuple(c[2:4]*1e-10),
                       exchanger_effectiveness=p0.exchanger_effectiveness if args.baseline else c[4])

    started = time.perf_counter()
    calls = 0
    best = np.inf

    def residual(x):
        nonlocal calls,best
        p = parameters(x[6*n:])
        rows = []
        for state,plate,load in zip(x[:6*n].reshape(6,n),plates,heater):
            balances,values,_ = ex.evaluate(state,plate,p,args.alpha)
            rows.append(np.r_[balances*1e4,values[:2]-load])
        r = np.concatenate(rows)
        calls += 1
        if calls % 500 == 1:
            print("progress",calls,round(time.perf_counter()-started,1),"norm",float(r@r),flush=True)
        return r

    cap = min(ex.drive.relief_swing)/1e6
    lower = np.r_[np.tile(np.r_[np.full(n-2,np.log(2)),cap/2+.01,.001],6),
                  [1e-7,1e-7,0,0], [] if args.baseline else [1e-5]]
    upper = np.r_[np.tile(np.r_[np.full(n-2,np.log(300)),2.5,cap],6),
                  [10,10,100,100], [] if args.baseline else [100]]
    initial = np.clip(np.r_[states.ravel(),coefficients],lower+1e-10,upper-1e-10)
    result = least_squares(residual,initial,bounds=(lower,upper),x_scale="jac",
        diff_step=1e-5,ftol=1e-9,xtol=1e-9,gtol=1e-9,max_nfev=args.max_nfev)
    p = parameters(result.x[6*n:])
    rows = []
    for i,state in enumerate(result.x[:6*n].reshape(6,n)):
        balance,values,details = ex.evaluate(state,plates[i],p,args.alpha)
        row = dict(point=i,approximate_balance_W=float(max(abs(details["heat"][ex.nodes]))))
        try:
            state,values,details = ex.solve(state,plates[i],p,args.alpha)
            row.update(state=state.tolist(),temperature_K=details["temperature"].tolist(),
                residual_W=(heater[i]-values[:2]).tolist(),diagnostics=values.tolist(),
                exact_balance_W=float(max(abs(details["heat"][ex.nodes]))))
            if i in (1,4):
                row["sensitivity"] = ex.sensitivity(state,plates[i],alpha=args.alpha,p=p)
        except RuntimeError as error:
            row["error"] = str(error)
        rows.append(row)
    output = dict(baseline=args.baseline,alpha=args.alpha,parameters=asdict(p),
        optimizer_success=result.success,optimizer_message=result.message,nfev=result.nfev,
        elapsed_s=time.perf_counter()-started,points=rows)
    if all("residual_W" in row for row in rows):
        output["rms_W"] = np.sqrt(np.mean(np.array([r["residual_W"] for r in rows])**2,axis=0)).tolist()
    suffix = "_stage2_seed" if args.seed_stage2 else ""
    path = HERE/f"quadrature_{'baseline' if args.baseline else 'current'}_{args.alpha:g}{suffix}.json"
    path.write_text(json.dumps(output,indent=2)+"\n")
    print("FINAL",json.dumps({k:v for k,v in output.items() if k!="points"}),flush=True)


if __name__ == "__main__":
    main()
