# T-World with intracellular ATP/ADP dependence: summary and lower-bound effect on the Ca²⁺ transient

*Branch `lowAtp`. Full write-up, figures and references: [README.md](README.md).*

## What was added to the model

All changes are opt-in. With the defaults, the model output is **bitwise identical** to
the current T-World for the endo, epi and mid cell types (checked with `isequal` on the
full trajectories).

- **New parameters**
  - `ATP_i` and free `ADP_i` (mM), clamped during a run.
  - The reference state `ATP_ref` = 8 mM and `ADP_ref` = 93.7 µM. `ADP_ref` is the human
    creatine kinase equilibrium (Weiss et al. 2005 PCr/ATP/total-Cr pools, K_CK = 1.66·10⁹
    M⁻¹, pH 7.05), computed by the new `getCKEquilibriumADP.m`.
  - `IKATP_Multiplier`, `KATP_Km_Multiplier` and `ca50_Multiplier`, all 1 by default.
- **Pump scaling (SERCA V_max, I_NaK, I_pCa).** The Cortassa et al. (2006) ATP/ADP factors,
  checked against the appendix (Eqs. A35, A46, A50).
  - Each factor is **divided by its value at the reference state**, so it is exactly 1
    there.
  - The I_pCa factor is a sum with maximum 2, so normalising it is essential.
- **I_KATP.** The human version of Ferrero et al. (1996) by Mena Tobar, Ferrero & Rodriguez
  Matas (CinC 2015): g₀ = 4.5, K_m = α(35.8 + 17.9·ADP^0.256) µM, H = 1.3 + 4.44·e^(−0.09·ADP),
  α = 0.32/0.86/1.0 for endo/mid/epi.
  - Only the open fraction **in excess of the reference** is added. The resting K_ATP
    conductance (~0.07 % open in endo) is taken to be absorbed in the calibrated background
    K⁺ currents.
  - I_KATP enters I_K_tot, so it acts on both V_m and [K⁺]_i.
- **`ca50_Multiplier`** scales the Land model TnC ca50, to probe myofilament Ca²⁺
  sensitivity.
- **New outputs.** `getCurrentsStructure` gains `IKATP` and the four factors `fATP_*`
  (output indices 36–40).
  - `model_TWorld` takes a trailing `atpParams` argument, which defaults when it is
    omitted.
  - Existing callers of the model are unaffected.

## Question and approach

What is the **smallest** change in the Ca²⁺ transient (CaT) when ATP falls from 8 mM to
4 or 2 mM, given current data? Every uncertain choice was set in the effect-reducing
direction:

- **Free ADP.** The lowest published ischaemic free ADP at that fraction of ATP: 182 µM
  at 50 % and 109 µM at 25 % (guinea pig; Ferrero et al. 2023 / Terkildsen et al. 2007).
  - Rat data (Bak & Ingwall 1994): 284 µM and 767 µM.
  - CK equilibrium cannot predict these values, because PCr is exhausted before ATP
    halves.
- **ADP_ref.** The human value of 94 µM, which is high compared with 15–39 µM in animals.
- **Cell type.** Endocardial, the least K_ATP-sensitive (α = 0.32).
- **ADP in the pump factors.** Free ADP, although Cortassa's constants use total ADP.
- **Omitted mechanisms.** Mechanisms that would only increase the effect are omitted:
  - Pi in SERCA thermodynamics;
  - MgADP myofilament sensitisation;
  - ATP activation of RyR2.

Protocol: endo cell at 1 Hz, nucleotides clamped, last of 300 beats (150 beats for
sweeps).

## Results

Changes are relative to 8 mM at matched myofilament sensitivity.

| | 4 mM ATP | 2 mM ATP |
|---|---|---|
| **CaT amplitude, lower bound** | **−5.8 %** | **−13.1 %** |
| CaT amplitude, I_KATP g₀ halved | −3.5 % | −7.0 % |
| CaT amplitude, I_KATP blocked (pumps only) | −1.1 % | −0.3 % |
| CaT amplitude, highest published ADP | −7.6 % | −25.5 % |
| APD₉₀ (lower bound) | −8.2 % | −21.0 % |
| Peak active tension (lower bound) | −6.9 % | −16.7 % |
| CaT decay τ (lower bound) | +0.4 % | −3.8 % |
| Excess K_ATP open fraction | 0.13 % | 0.36 % |
| f_SERCA / f_NaK / f_pCa | 0.980 / 0.996 / 0.986 | 0.990 / 0.994 / 0.959 |

1. **The effect is almost entirely electrical.**
   - Opening a fraction of a percent of K_ATP channels shortens the AP. That lowers L-type
     Ca²⁺ entry and, slightly, SR load (−1 %).
   - The pumps remain >96 % saturated, so the direct effect of ATP on Ca²⁺ handling is
     ≈1 %.
2. **Decay kinetics are not a robust marker.** In the lower-bound scenarios, τ, CaTD₅₀ and CaTD₉₀ change by less than
   4 %; even at the highest published ADP, CaTD₉₀ changes by at most +7 %.
3. **Myofilament Ca²⁺ sensitivity is independent of the ATP effect but can mask it.**
   - Across ca50 ×0.9–2, the relative ATP effect stays at −5.8 to −6.1 % (4 mM) and −13.1
     to −13.7 % (2 mM).
   - The absolute CaT is strongly ca50-dependent through TnC buffering: ca50 ×1.25 gives
     +14 % amplitude and −13 % τ.
   - Desensitisation of about ×1.1 (4 mM) or ×1.27 (2 mM) is enough to cancel the
     absolute reduction. Acidosis desensitises more than this (Palmer & Kentish 1994).
     The observable CaT in acidotic ischaemia may therefore not fall, even though the
     ATP-attributable component does.

## Limitations

- **The bound depends on K_ATP.**
  - Without functional K_ATP channels, the guaranteed effect is ≈1 %.
  - The human gating rests on recombinant channel data: Babenko et al. (1998), IC₅₀ ≈
    20 µM and Hill ≈ 1. Our H = 1.3 is the conservative choice.
  - The parameters come from a 4-page proceedings paper, fitted in TP06. g₀ is printed
    as "mS"; mS/µF is inferred.
  - The Mg²⁺/Na⁺ block and Q₁₀ constants of Ferrero et al. (1996) are from secondary
    sources and are unverified.
  - The rise in free Mg²⁺ as MgATP falls is not modelled. It would strengthen the Mg²⁺
    block, which is partly covered by the g₀ ×0.5 row.
  - Sanity check: Ferrero et al. (1996) report ~0.6 % opening → 50 % APD shortening.
    Here 0.36 % gives 21 %, so the current is not overstated.
- **ADP data.** All ischaemic ADP values are animal and CK-derived, and they are
  semi-quantitative at acidic pH. Cortassa's constants are in-vitro, guinea pig/rat.
- **Not modelled:**
  - Pi;
  - extracellular K⁺ accumulation;
  - direct acidosis effects on channels and exchangers;
  - L-type rundown;
  - time-varying ischaemia (states are clamped);
  - [Mg²⁺]_i changes.
- **[Na⁺]_i drift.** [Na⁺]_i still drifts by 0.1 mM at 300 beats; the CaT is stable after
  about 50 beats.
- **Falsifiable prediction.**
  - In human myocytes dialysed with 2 mM vs 8 mM MgATP (~0.1 mM ADP, pH and Mg²⁺
    clamped), the CaT amplitude should fall by at least ~7–13 %.
  - Glibenclamide should abolish that fall.

## Questions for you (T-World)

- Is normalising to the reference state the right way to keep the calibrated baseline?
  The alternative is re-fitting the maximal rates under explicit nucleotide factors.
- Should the TP06-derived K_ATP g₀ be re-scaled to T-World's repolarisation reserve? If
  so, is there human APD-under-pinacidil (or metabolic-inhibition) data you would
  calibrate against?
