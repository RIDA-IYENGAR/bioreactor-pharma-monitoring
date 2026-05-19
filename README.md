# Bioreactor Monitoring and Control for Pharmaceutical Production

**MathWorks MATLAB & Simulink Challenge Project**
Applying Quality by Design (QbD) and Process Analytical Technology (PAT)
to industrial-scale penicillin fermentation using IndPenSim.

---

## Team

| Name | GitHub | Contribution |
|------|--------|--------------|
| Rida Iyengar | [@RIDA-IYENGAR](https://github.com/RIDA-IYENGAR) | Phases 1, 2, 3 — Simulation, Multivariate Analysis, Soft Sensor |
| Saniha M Shetty | [@sanihamshetty](https://github.com/sanihamshetty) | Phase 4 — pH PID Controller |
| Nandana Biju | [@Nandanabiju24btech](https://github.com/Nandanabiju24btech) | Phase 5 — Substrate Feed Control and Integration |

---

## Project Overview

Biopharmaceutical manufacturing suffers from batch-to-batch variability in
penicillin yield, causing quality inconsistencies and financial losses at
industrial scale. This project builds a full QbD/PAT pipeline on the IndPenSim
industrial penicillin fermentation simulator to:

1. Quantify natural batch variation across multiple simulated runs
2. Identify which process variables drive yield inconsistency (CPPs)
3. Predict penicillin concentration in real-time using a Raman spectroscopy soft sensor
4. Control the top CPPs using PID feedback controllers

---

## Repository Structure
bioreactor-pharma-monitoring/
├── code/
│   ├── first_run.m                    # Phase 1: Baseline IndPenSim run
│   ├── multi_batch_run.m              # Phase 1: 20-batch variation study
│   ├── phase2_CPP_analysis.m          # Phase 2: PCA and PLS analysis
│   ├── phase3_raman_batch.m           # Phase 3: Raman data collection
│   ├── phase3_soft_sensor.m           # Phase 3: Soft sensor training
│   └── [IndPenSim simulator files]
├── phase4_pH_control/
│   └── phase4_pH_PID_controller.m     # Phase 4: pH PID controller
├── phase5_substrate_control/
│   └── phase5_substrate_PID_controller.m  # Phase 5: Dual-loop controller
├── results/                           # All generated plots and data files
├── docs/
└── README.md

---

## Requirements

- MATLAB R2021a or later
- Statistics and Machine Learning Toolbox
- Control System Toolbox
- Simulink
- IndPenSim V2 (included in code/ folder)

---

## How to Run

```matlab
% Step 1: Open MATLAB and navigate to the code/ folder

% Step 2: Add all files to path
addpath(genpath(pwd))
addpath(genpath('../phase4_pH_control'))
addpath(genpath('../phase5_substrate_control'))

% Step 3: Run phases in order (wait for each to finish)
run('first_run.m')                  % ~2 min
run('multi_batch_run.m')            % ~20 min
run('phase2_CPP_analysis.m')        % ~1 min
run('phase3_raman_batch.m')         % ~15 min
run('phase3_soft_sensor.m')         % ~5 min
run('../phase4_pH_control/phase4_pH_PID_controller.m')            % ~5 min
run('../phase5_substrate_control/phase5_substrate_PID_controller.m')  % ~10 min
```

All plots and results are saved automatically to the `results/` folder.

---

## Results

### Phase 1 — Batch Variation Study (20 batches)

| Metric | Value |
|--------|-------|
| Mean penicillin yield | 24.22 g/L |
| Standard deviation | 7.38 g/L |
| Coefficient of variation | 30% |
| Min yield | 13.67 g/L |
| Max yield | 34.45 g/L |
| Biomass SD (for comparison) | 1.06 g/L |

Penicillin varies 7x more than biomass, confirming it as the
Critical Quality Attribute (CQA) requiring active control.

---

### Phase 2 — CPP Identification via PCA and PLS

| Rank | Critical Process Parameter | PLS Coefficient |
|------|---------------------------|-----------------|
| 1 | Substrate variability (SD) | 0.517 |
| 2 | pH variability (SD) | 0.329 |
| 3 | Substrate mean level | 0.288 |
| 4 | pH mean level | 0.203 |
| 5 | Dissolved oxygen mean | 0.177 |

PLS model R² = 0.944 — 94.4% of penicillin yield variation is
predictable from process measurements alone.
PCA: PC1 + PC2 explain 70.9% of total process variance.

---

### Phase 3 — Raman Spectroscopy Soft Sensor

| Metric | Value |
|--------|-------|
| Training data size | 10 batches × 1150 timepoints × 2200 wavelengths |
| PLS components | 15 |
| R² on test batches | 0.9998 |
| RMSEP | 0.1354 g/L |
| Prediction frequency | Every 12 minutes |
| Off-line assay delay replaced | 4+ hours |

The soft sensor reduces measurement latency from 4+ hours to 12 minutes
with 0.4% relative prediction error, enabling true real-time monitoring.

---

### Phase 4 — pH PID Controller
*Contributor: Saniha M Shetty*

PID feedback controller maintaining pH at setpoint (6.5) throughout
the batch, targeting the second-ranked CPP from Phase 2.

---

### Phase 5 — Substrate Control and Dual Loop Integration
*Contributor: Nandana Biju*

PID controller for substrate feed rate (top CPP) combined with the
Phase 4 pH controller into a full dual-loop integrated control system.

---

## Methodology Flow
IndPenSim Simulator (industrial penicillin fermentation)
|
v
Phase 1: 20-batch study → quantify natural variation (SD = 7.38 g/L)
|
v
Phase 2: PCA + PLS → identify CPPs (Substrate SD, pH SD)
|
v
Phase 3: Raman spectra → PLS soft sensor (R² = 0.9998)
|
v
Phase 4: pH PID controller → reduce pH variability
|
v
Phase 5: Substrate PID + Dual-loop → full process control

---

## References

1. Goldrick et al. (2015). The development of an industrial-scale
   fed-batch fermentation simulation. Journal of Biotechnology.
2. Goldrick et al. (2019). Modern monitoring and control of industrial
   fermentation processes. Computers and Chemical Engineering.
3. MathWorks IndPenSim V2 — MATLAB File Exchange (ID: 49041)

---

## Institution

PES University, Bengaluru — Department of Biotechnology Engineering
MathWorks MATLAB and Simulink Challenge 2024-2025
