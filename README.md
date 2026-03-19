# Bioreactor Monitoring and Control for Pharmaceutical Production

**MathWorks MATLAB & Simulink Challenge Project**  
**Student:** [Your Name] | PES University, Biotechnology  

## Project Overview
Implementing Quality by Design (QbD) and Process Analytical Technology (PAT) 
on an industrial penicillin fermentation process using IndPenSim in MATLAB/Simulink.

## Progress
- [X] Phase 1: IndPenSim familiarization and first run
- [X] Phase 2: CPP/CQA identification via multivariate analysis
- [ ] Phase 3: Soft sensor development using Raman spectroscopy
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
## Key Results So Far
**Phase 1:** Natural batch variation SD = 7.38 g/L (CV = 30%) across 20 batches.

**Phase 2 CPPs identified (PLS R² = 0.944):**
1. Substrate variability (SD) — top driver of yield variation
2. pH variability (SD) — second most critical
3. Substrate mean level — third most critical

**CQA:** Final penicillin concentration (g/L)
