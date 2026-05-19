%% =========================================================
% run_full_project.m
% MASTER SCRIPT — Runs all phases in sequence
% Bioreactor Monitoring and Control — MathWorks Challenge
%
% Team:
%   Rida Iyengar    — Phases 1, 2, 3
%   Saniha M Shetty — Phase 4
%   Nandana Biju    — Phase 5
%
% Run this script to reproduce all results from scratch.
% Expected total runtime: 45-60 minutes
% =========================================================
clc; clear; close all;

% Add all project files to path
addpath(genpath(pwd))
addpath(genpath('../phase4_pH_control'))
addpath(genpath('../phase5_substrate_control'))

fprintf('==============================================\n')
fprintf(' BIOREACTOR MONITORING & CONTROL PROJECT\n')
fprintf(' MathWorks MATLAB & Simulink Challenge\n')
fprintf('==============================================\n\n')

%% --- PHASE 1: Baseline simulation ---
fprintf('>>> PHASE 1: Running baseline batch...\n')
run('first_run.m')
fprintf('>>> PHASE 1 COMPLETE\n\n')

%% --- PHASE 1 Step 2: Multi-batch variation ---
fprintf('>>> PHASE 1 Step 2: Running 20 batches...\n')
run('multi_batch_run.m')
fprintf('>>> PHASE 1 Step 2 COMPLETE\n\n')

%% --- PHASE 2: CPP/CQA Analysis ---
fprintf('>>> PHASE 2: PCA and PLS analysis...\n')
run('phase2_CPP_analysis.m')
fprintf('>>> PHASE 2 COMPLETE\n\n')

%% --- PHASE 3: Raman data collection ---
fprintf('>>> PHASE 3 Step 1: Collecting Raman batch data...\n')
run('phase3_raman_batch.m')
fprintf('>>> PHASE 3 Step 1 COMPLETE\n\n')

%% --- PHASE 3: Soft sensor training ---
fprintf('>>> PHASE 3 Step 2: Training soft sensor...\n')
run('phase3_soft_sensor.m')
fprintf('>>> PHASE 3 Step 2 COMPLETE\n\n')

%% --- PHASE 4: pH PID Controller ---
fprintf('>>> PHASE 4: pH PID Controller...\n')
run('../phase4_pH_control/phase4_pH_PID_controller.m')
fprintf('>>> PHASE 4 COMPLETE\n\n')

%% --- PHASE 5: Substrate + Dual Loop Controller ---
fprintf('>>> PHASE 5: Substrate PID + Dual Loop Controller...\n')
run('../phase5_substrate_control/phase5_substrate_PID_controller.m')
fprintf('>>> PHASE 5 COMPLETE\n\n')

%% --- Final Summary ---
fprintf('==============================================\n')
fprintf(' ALL PHASES COMPLETE\n')
fprintf('==============================================\n')
fprintf('\n--- RESULTS SUMMARY ---\n')
fprintf('Phase 1 : Baseline batch complete\n')
fprintf('          20-batch SD = 7.38 g/L (CV=30%%)\n')
fprintf('Phase 2 : Top CPPs = Substrate SD, pH SD\n')
fprintf('          PLS R2 = 0.944\n')
fprintf('Phase 3 : Soft sensor R2 = 0.9998\n')
fprintf('          RMSEP = 0.1354 g/L\n')
fprintf('Phase 4 : pH PID controller — see results/\n')
fprintf('Phase 5 : Dual-loop controller — see results/\n')
fprintf('\nAll plots saved to results/ folder\n')
fprintf('==============================================\n')
