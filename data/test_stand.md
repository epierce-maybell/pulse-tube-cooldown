# Pulse-tube test-stand records

The test stand runs one pulse-tube head at a time and logs three Lake Shore Cernox thermometers: T1 on stage 1, TPR on the regenerator clamp, and T2 on stage 2. Each `<system>/pulse_tube/` folder holds that head's cooldown traces (`cooldown_<date>.dat`) and load maps (`load_map_<date>.dat`). The stand's full archive, 43 continuous logs and 22 load maps across eleven heads from March to October 2026, is in `../../data/test_stand/` relative to the repository.

## Format

A cooldown file opens with the calibration curve for each channel, T2 through Z1512 and T1 and TPR both through Z1518, then a column line. Each row holds the timestamp, elapsed minutes, T1, TPR, T2, two copies of T1 (`Taux` and `Taux2`), and unused fields. Rows come every 2 s, and TPR updates about every 15 s. The elapsed column can restart within a file, so the timestamp is the clock. `sn09/pulse_tube/cooldown_2026-05-14.dat` has a second short run appended from 14:19:17.

A load map opens with the valve settings, the valve frequency (`f:1.1`), and the applied powers, then holds one row per steady state: stage-1 and stage-2 heater power in W, elapsed minutes, T1, TPR, T2, `Taux`, and `Taux2`.

## Channel behavior

Found across twelve cooldowns of at least six heads, all logged with the same curves.

- **T1** reads NaN between 70 and 76 K and between 219 and 230 K. Above 230 K it latches at a stale value, 259 to 292 K at rest.
- **T2** reads NaN between 70 and 76 K, and between 20 and 21 K on some heads, SN09 included. At rest it reads 306 to 310 K while TPR reads 292 to 296 K.
- **TPR** reads 33 to 35 K at base on every head, and SN09's two load maps agree within 0.7 K at matching loads. During cooldown it rises in a pattern that repeats on a given head and differs between heads: above room temperature on SN09 (325 and 351 K), SN038 (365 and 376 K), and SN035 (385 K), and below it on SN007, SN019, and SN028. It then drops to base within about 35 min once T1 reaches 45 to 60 K. The curves are identical across heads, so the rise is not a calibration error. Its cause is unknown. It is repeatable but does not read the clamp: the clamp sits 0.658 of the way from stage 1 to stage 2, yet across SN09's load-map points only 8% to 38% of the T1-to-T2 drop happens above TPR, while on the SN09 system the RGP plate clamped at the same position reads 11.1 K at base with PRP at 72.6 K and CFP at 2.8 K, 88% of the way down.

## Use

Load-map fits use T1 and T2 at each point. Cooldown comparisons use T1 below 219 K and T2 throughout. TPR is a diagnostic to plot beside the model's clamp node, not a target; the clamp is tested by the SN09 system's RGP trace.
