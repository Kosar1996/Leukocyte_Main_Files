# Reflection 10/14: code changes since the Sep 8 handover

**Baseline:** previous version = `Leukocyte_Main_Files_Backup` @ `d692f23` (Sep 8 handover; repository history from `683c85d` "Restore Sep 11 code version", Sep 14, to `d692f23`, Sep 17).
**Updated code:** `Leukocyte_Main_Files` (`524179a` → `4a4ac3e`) and `Leukocyte-1002-Serial` (to `eb8d3d6`, `ad86de3`) / `-Parallel` (to `8df6d1e`, `20e01e9`), 10/13.
**Classes:** **N** = numerical/calculation mechanics · **P** = physical equation / continuum-mechanics formulation · **S** = model setup.

## 1–2. Changes compared with previous version, classified

| Package (commit; Serial repo from 10/1) | Change | Class |
|---|---|---|
| 9/21 (`524179a`) | Endothelium assembler: neo-Hookean deviatoric stress `J^(-5/3)` → `J^(-2/3)` | P |
| 9/21 (`524179a`) | Interface wall velocity: old displacement evaluated at the new interface position (material velocity) | P |
| 9/21 (`524179a`) | Standalone solid solve (used by the traction correction): separate Kelvin–Voigt call removed and `uOld` passed to the assembler instead; the assembler takes 3 inputs, so the call fails with "Too many input arguments" and the traction correction is skipped | P |
| 9/21 (`524179a`) | Curve interpolation clamped at the ends; Newton safeguards in fluid-stress point location; stress recovery at 4 Gauss points | N |
| 9/22 (`524179a`) | Solid Newton: 5 → 1 load step, LM shift instead of fsolve fallback, best iterate on failure, new tolerances; `uOld` passed via `par` but not read by the assembler | N |
| 9/23 (`b435738`) | Outer radius 4 → 6 µm; leukocyte prestress P600 → P300; full 2D body-fitted MAC fluid instead of the hybrid gap model; global 2D pressure traction | S |
| 9/23 (`b435738`) | Soft pressure saturation (tanh) above 50 kPa | N |
| 9/28 (`60f34f3`) | Traction correction relaxation 1.0 → 0.05, tol 1e-4 → 1e-2; new outer Aitken traction-coupling loop (ω 0.08–0.50); dt 1e-5 → 2.5e-6 | N |
| 9/28 (`60f34f3`) | Limiters: 5 kPa traction cap, 10 mm/s wall-velocity cap, non-converged solid keeps 10 % of step, monolithic 50 % | N |
| 9/28 (`60f34f3`) | Interface velocity back to the Eulerian difference; flux-equivalent radial wall velocity | P |
| 9/29 (`8d6e510`) | Viscous stress `2·η·D` (η = 10 Pa·s) added in the solid assembler, without its tangent term | P |
| 9/29 (`8d6e510`) | Monolithic: LM fallback, Armijo line search, error instead of partial acceptance; caps p ≤ 3000 Pa, shear ≤ 500 Pa | N |
| 10/1 (`9e31ea1`) | Active translocation body force `f0` on the leukocyte, tanh gate at the pore | P |
| 10/1 (`9e31ea1`) | Error if J_min ≤ 0.05; line search must reach 0.8·R0; gap floor `δE ≥ δL + gapFloor` | N |
| 1002/1005 (`7ac67eb`, `9e31ea1`) | (mine) restart-state fix, speed-ups (same changes in the parallel version) | N |
| 1008 (`70ecc0a`) | Force gate uses min z for f0 < 0 | P |
| 1009 (`66b66fa`) | Raised-cosine force ramp (t_ramp = 3 µs); (mine) ramp starts at run start | N |
| 1013 (`eb8d3d6`) | (mine, merging the 10/8 stabilization files) dt cap at step 14, dt-halving retry, v_wall defaults; step 14 bit-identical when restarted from step 13 | N |

## 3. Reflection: structure, coupling, solver

**How the updates are designed to address the instabilities.**
- The under-relaxed traction correction and the outer Aitken loop damp traction oscillations between fluid and solid in the thin gap.
- A smaller dt, the step-14 cap and the dt-halving retry aim to keep each step within Newton's convergence range.
- When fsolve stalls, LM and the Armijo line search take over.
- The gap floor prevents interpenetration; the J_min check stops a step with an error when elements become nearly inverted.
- The force ramp avoids an impulsive start of `f0`.
- Caps, saturation and partial acceptance are intended to limit pressure spikes and prevent blow-up; they also mean a step can be accepted with the coupled equations satisfied only approximately.

**Continuum mechanics.** In the previous version, `J^(-2/3)` was already used for the leukocyte, the P300 initialisation and stress recovery. The 9/21 package applied it to the endothelium assembler but not to `assemble_finite_def_internal_force_only`, which still uses `J^(-5/3)`. The draft paper (eq. 19) writes the endothelium Cauchy stress as `G·J^(-5/3)·(B − tr(B)/3·I) + K(J−1)·I`; since the assembler uses `T` as Cauchy stress (`P = J·T·F^-T`), `J^(-2/3)` corresponds to the Kirchhoff form. The leukocyte in the draft is fractional Kelvin–Voigt with small strain (eq. 24); the code uses finite-strain neo-Hookean plus a Kelvin–Voigt term, as in the previous version.

Viscosity: the monolithic solve keeps Kelvin–Voigt with its tangent. The standalone solid solve had none from 9/21 to 9/28, and since 9/29 it uses a different viscous law without its tangent term. This may contribute to the stalled solid solves seen on 10/13 (1,882 of 1,968).

**Coupling and solver: open points.**
1. The physical-residual norms checked after fsolve (one of the two acceptance criteria) come from the unscaled residual. It uses `assemble_finite_def_internal_force_only` (`J^(-5/3)`) and applies `f0` without the tanh pore gate, unlike the scaled residual that fsolve solves.
2. In the retry block, the flag set is `acceptStep`, but the loop tests `acceptedStep`. After `break`, execution continues after the coupling loop instead of restarting the step.

**Next steps (proposed):**
1. One stress law in all code paths.
2. One viscous law, with its tangent term, in both solves.
3. The same pore gate in both residuals.
4. Fix the retry flag.
5. Then remove limiters one at a time against the step-28 reference.
