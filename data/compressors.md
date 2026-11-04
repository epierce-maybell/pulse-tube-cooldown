# Compressor records

Every compressor is a Sumitomo F-70H. The test stand runs one compressor on its head; SN09 and FP1 are in `systems.md`. A head's valve frequency is its shipped setting, a per-unit number rather than a drive value (see `test_stand.md`).

## Charge and pressures

- Every head in `docs/pt_performance_overview_2026-10-30.pdf` is filled to 17.0 to 17.1 bar at rest.
- The F-70H/F-70L/F-70LP Operating Manual (SHI 267472A, rev. E, 2011), for a 4K pulse-tube head: equalization pressure at 20 C of 13.6, 14.3, or 16.6 bar gauge nominal, without saying which configuration takes which; factory charge of 1430 to 1460 kPa at 20 C; operating supply pressure of 18 to 25 bar gauge, varying with head and load; shutdown on return pressure above 140 psig or below 15 psig; power draw of 6.7 to 7.2 kW on 50 Hz mains and 8.0 to 8.5 kW on 60 Hz. Neither that manual nor the F-70H/F-70L Service Manual gives a displacement rate or internal-relief setpoint.
- The running mean pressure falls as gas moves into the cold volumes. S0023's Münster cooldown trace (overview, p. 33) shows the return pressure rising from about 6.8 bar to 8.7 bar an hour in and settling near 6.5 bar at base.
- Front panel of one of FP1's compressors during cooldown, with the head at 260 to 280 K and the swing at its relief cap: supply 360 to 380 psig on the analog gauge, return 90 psig on the LCD. That puts the relief swing at 270 to 290 psi (18.6 to 20.0 bar) and the running mean at 230 psig (16.9 bar absolute), closest to the manual's 16.6 bar gauge equalization option.

## Displacement rate

The scroll capsule is a Copeland ZCH72C3G/C4G, rated 3500 rpm on 60 Hz. Copeland does not publish its displacement; the ZR72KC on what appears to be the same frame sweeps 5.98 in^3/rev, or 20.6 m^3/h (5.72e-3 m^3/s) at 3500 rpm. The swept rate bounds the displacement rate from above, since helium leakage, warm suction gas, and injected oil reduce delivery. It stands as the displacement rate until a reading at base replaces it, and delivers 6.6 g/s at the 90 psig return and 300 K.

## Manifold volume

On the test stand, the warm volume is two 20 m DN25 flexlines and two 5 L buffers. That is the default for the fridges until measured.
