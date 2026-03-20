# Bioreactor Monitoring and Control for Pharmaceutical Production

**MathWorks MATLAB & Simulink Challenge Project**  
## Team
| Name | Role |
|------|------|
| Rida Iyengar | Phase 1, 2, 3 — Simulation, Multivariate Analysis, Soft Sensor |
| Saniha M Shetty | Phase 4 — pH PID Controller |
| Nandana Biju | Phase 5 — Substrate Feed Control & Integration |

## Project Overview
Implementing Quality by Design (QbD) and Process Analytical Technology (PAT) 
on an industrial penicillin fermentation process using IndPenSim in MATLAB/Simulink.

## Progress
- [X] Phase 1: IndPenSim familiarization and first run
- [X] Phase 2: CPP/CQA identification via multivariate analysis
- [X] Phase 3: Soft sensor development using Raman spectroscopy
- [ ] Phase 4: PID control design in Simulink

## How to Run
Instructions will be added as the project progresses.

## References
- Goldrick et al. (2014) - IndPenSim simulator development
- Goldrick et al. (2019) - Monitoring and control of IndPenSim

## Key Results So Far
**Phase 1 Finding:** Natural batch-to-batch variation in penicillin yield 
shows SD = 7.38 g/L (CV = 30%) across 20 batches, with yield ranging from 
13.67 to 34.45 g/L. Biomass variation is comparatively low (SD = 1.06 g/L), 
confirming penicillin production is the critical quality attribute requiring control.

**Phase 2 CPPs identified (PLS R² = 0.944):**
1. Substrate variability (SD) — top driver of yield variation
2. pH variability (SD) — second most critical
3. Substrate mean level — third most critical

**CQA:** Final penicillin concentration (g/L)

**Phase 3 — Soft Sensor (Raman PLS Model):**
- 10 batches × 1150 timepoints × 2200 wavelengths = 11,500 training samples
- PLS model with 15 components
- R² = 0.9998 on held-out test batches (Batches 8–10)
- RMSEP = 0.1354 g/L (vs concentration range 0–35 g/L)
- Enables real-time penicillin prediction every 12 minutes
  vs 4+ hour off-line lab delay

| Phase 4 | In progress |
| Phase 5 | In progress |
