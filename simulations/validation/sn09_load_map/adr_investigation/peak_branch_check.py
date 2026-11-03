"""Try both sides of the GOS heat-capacity peak at a failed corrected-drive fit."""

import json

import numpy as np
from scipy.optimize import least_squares

from flow_consistency import Experiment, HERE, fit


def main():
    source = HERE/"continuation_baseline_drive_loss_seed.json"
    saved = json.loads(source.read_text())
    # An unsuccessful initial refit has no verified final state.
    saved["equilibrium_verified"] = False
    saved["points"] = []
    source.write_text(json.dumps(saved,indent=2)+"\n")
    attempt = saved["history"][-1]
    ex = Experiment(True,True)
    p = fit.Parameters(**attempt["parameters"])
    state = np.array(attempt["failed_states"][0])
    plate = np.array(ex.saved["points"][0]["stage_temperature_K"])
    rows = []
    for t in (4.8,5.2,5.24,5.26,5.4,6,8,12):
        guess = state.copy()
        guess[-3] = np.log(t)
        try:
            root,values,details = ex.solve(guess,plate,p,1)
            row = dict(initial_GOS_K=t,final_GOS_K=float(np.exp(root[-3])),
                cooling_W=values[:2].tolist(),state=root.tolist(),
                max_balance_W=float(max(abs(details["heat"][ex.nodes]))))
        except RuntimeError as error:
            row = dict(initial_GOS_K=t,error=str(error))
        rows.append(row)
        print(json.dumps(row),flush=True)
    (HERE/"peak_branch_check.json").write_text(json.dumps(rows,indent=2)+"\n")
    root = np.array(next(r["state"] for r in rows if "state" in r))
    free = np.delete(np.arange(len(root)),len(root)-3)
    cap = min(ex.drive.relief_swing)/1e6
    lower = np.r_[np.full(len(root)-2,np.log(2)),cap/2+.01,.001][free]
    upper = np.r_[np.full(len(root)-2,np.log(300)),2.5,cap][free]
    scan = []
    for t in (10,8.5,7.81,7,6,5.4,5.26,5.25,5.245,5.24,5.2,5,4.8,4.3,4):
        root[-3] = np.log(t)

        def residual(y):
            candidate = root.copy()
            candidate[free] = y
            return ex.evaluate(candidate,plate,p,1)[0][free]

        solved = least_squares(residual,root[free],bounds=(lower,upper),jac="3-point",
            x_scale="jac",ftol=1e-12,xtol=1e-12,gtol=1e-12,max_nfev=150)
        root[free] = solved.x
        residuals,values,details = ex.evaluate(root,plate,p,1)
        row = dict(GOS_K=t,GOS_net_heat_W=float(details["heat"][-2]),
            other_equation_error=float(max(abs(residuals[free]))),cooling_W=values[:2].tolist())
        scan.append(row)
        print("SCAN",json.dumps(row),flush=True)
    (HERE/"peak_balance_scan.json").write_text(json.dumps(scan,indent=2)+"\n")


if __name__ == "__main__":
    main()
