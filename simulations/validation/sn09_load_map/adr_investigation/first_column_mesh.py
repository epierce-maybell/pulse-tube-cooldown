"""Refine stage 1 at fixed stage-2 fit ratios, refitting linear coefficients."""

from dataclasses import asdict, replace
import json

import numpy as np

from flow_consistency import Experiment, HERE, fit


def main():
    ex = Experiment()
    base, p = ex.head, ex.parameters
    rows = []
    for count in (1,2,4,8):
        head = replace(base,
            material=(tuple(m for m in base.material[0] for _ in range(count)),base.material[1]),
            mass=(np.repeat(base.mass[0]/count,count),base.mass[1]),
            pore_volume=(np.repeat(base.pore_volume[0]/count,count),base.pore_volume[1]),
            clamp=base.clamp+len(base.mass[0])*(count-1))
        old_x = np.r_[0,np.arange(len(base.mass[0]))+.5,len(base.mass[0])]
        new_x = (np.arange(len(head.mass[0]))+.5)/count
        guesses = []
        for saved in ex.saved["points"]:
            t = np.array(saved["temperature_K"])
            guesses.append(np.r_[t[0],np.interp(new_x,old_x,t[:base.stage_1+1]),t[base.stage_1:]])
        try:
            result = fit.fit_at_ratio(head,ex.drive,ex.walls,ex.data,warm_temperature=ex.warm,
                ratio=p.conductance[1]/p.packing_constant[1],
                effectiveness_ratio=p.exchanger_effectiveness/p.packing_constant[1],
                frequency=p.frequency,guesses=guesses)
        except RuntimeError as error:
            rows.append(dict(nodes_per_original_layer=count,error=str(error)))
            break
        ex.head,ex.nodes = head,fit.column_indices(head)
        slopes = []
        for i in (1,4):
            c = result.columns[i]
            state = np.r_[np.log(c.temperature[ex.nodes]),
                (c.supply_pressure+c.return_pressure)/2e6,(c.supply_pressure-c.return_pressure)/1e6]
            slope = ex.sensitivity(state,np.array(ex.saved["points"][i]["stage_temperature_K"]),p=result.parameters)
            slopes.append(slope["estimates"][-1]["dT1_dT2_at_fixed_Q1"])
        row = dict(nodes_per_original_layer=count,parameters=asdict(result.parameters),
            rms_W=np.sqrt(np.mean(result.residual**2,axis=0)).tolist(),cross_slopes=slopes,
            maximum_balance_W=max(float(max(abs(c.heat_flow[ex.nodes]))) for c in result.columns))
        rows.append(row)
        print(json.dumps(row),flush=True)
    (HERE/"first_column_mesh.json").write_text(json.dumps(rows,indent=2)+"\n")


if __name__ == "__main__":
    main()
