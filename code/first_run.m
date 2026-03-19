%% =========================================================
% first_run.m
% Phase 1: First IndPenSim run - Baseline batch
% =========================================================
clc; clear; close all;

%% --- Simulation Settings ---
h = 0.2;              % sampling time (hours) = 12 mins
T = 230;              % batch duration (hours) = standard penicillin batch
Batch_no = 1;         % this is batch number 1
Batch_time = 0:h:T;   % time vector

%% --- Control Flags ---
Ctrl_flags.SBC        = 0;  % standard batch control
Ctrl_flags.PRBS       = 0;  % recipe driven (no operator control)
Ctrl_flags.Fixed_Batch_length = 0;  % fixed batch length
Ctrl_flags.IC         = 0;  % randomised initial conditions
Ctrl_flags.Inhib      = 2;  % full inhibition (DO2, T, pH, CO2, PAA, N)
Ctrl_flags.Dis        = 1;  % process disturbances ON
Ctrl_flags.Faults     = 0;  % no faults for baseline
Ctrl_flags.Vis        = 0;  % simulated viscosity
Ctrl_flags.Raman_spec = 0;  % no Raman spectral data yet
Ctrl_flags.Batch_Num  = Batch_no;
Ctrl_flags.Off_line_m     = 12;  % off-line measurement every 12 hrs
Ctrl_flags.Off_line_delay =  4;  % 4 hr analysis delay
Ctrl_flags.plots      = 0;  % we'll make our own plots
Ctrl_flags.T_sp       = 298;   % temperature setpoint (K) = 25 C
Ctrl_flags.pH_sp      = 6.5;   % pH setpoint

%% --- Random Seed (for reproducibility) ---
Random_seed_ref = 42;   % fixed seed so results are same every run
Seed_ref = 31 + Random_seed_ref;
Rand_ref = 1;

%% --- Initial Conditions ---
rng(Seed_ref + Batch_no + Rand_ref); Rand_ref = Rand_ref+1;
intial_conds = 0.5 + 0.05*randn;

rng(Seed_ref + Batch_no + Rand_ref); Rand_ref = Rand_ref+1;
x0.mux = 0.41 + 0.025*randn;       % max specific growth rate of biomass [h^-1]

rng(Seed_ref + Batch_no + Rand_ref); Rand_ref = Rand_ref+1;
x0.mup = 0.041 + 0.0025*randn;     % max specific growth rate of penicillin [h^-1]

rng(Seed_ref + Batch_no + Rand_ref); Rand_ref = Rand_ref+1;
x0.S   = 1 + 0.1*randn;            % initial substrate concentration [g/L]

rng(Seed_ref + Batch_no + Rand_ref); Rand_ref = Rand_ref+1;
x0.DO2 = 15 + 0.5*randn;           % initial dissolved oxygen [mg/L]

rng(Seed_ref + Batch_no + Rand_ref); Rand_ref = Rand_ref+1;
x0.X   = intial_conds + 0.1*randn; % initial biomass [g/L]
x0.P   = 0;                         % initial penicillin [g/L]

rng(Seed_ref + Batch_no + Rand_ref); Rand_ref = Rand_ref+1;
x0.V   = 5.800e+04 + 500*randn;    % initial volume [L]

rng(Seed_ref + Batch_no + Rand_ref); Rand_ref = Rand_ref+1;
x0.Wt  = 6.2e+04 + 500*randn;      % initial weight [kg]

rng(Seed_ref + Batch_no + Rand_ref); Rand_ref = Rand_ref+1;
x0.CO2outgas = 0.038 + 0.001*randn; % CO2 offgas [%]

rng(Seed_ref + Batch_no + Rand_ref); Rand_ref = Rand_ref+1;
x0.O2  = 0.20 + 0.05*randn;        % O2 offgas [%]

rng(Seed_ref + Batch_no + Rand_ref); Rand_ref = Rand_ref+1;
x0.pH  = 6.5 + 0.1*randn;          % initial pH

rng(Seed_ref + Batch_no + Rand_ref); Rand_ref = Rand_ref+1;
x0.T   = 297 + 0.5*randn;          % initial temperature [K]

x0.a0  = intial_conds*(1/3);        % type a0 biomass [g/L]
x0.a1  = intial_conds*(2/3);        % type a1 biomass [g/L]
x0.a3  = 0;
x0.a4  = 0;
x0.Culture_age = 0;                 % culture age [hr]

rng(Seed_ref + Batch_no + Rand_ref); Rand_ref = Rand_ref+1;
x0.PAA = 1400 + 50*randn;           % initial PAA concentration [mg/L]

rng(Seed_ref + Batch_no + Rand_ref); Rand_ref = Rand_ref+1;
x0.NH3 = 1700 + 50*randn;           % initial nitrogen [mg/L]

rng(Seed_ref + Batch_no + Rand_ref); Rand_ref = Rand_ref+1;
alpha_kla = 85 + 10*randn;          % kla constant

rng(Seed_ref + Batch_no + Rand_ref); Rand_ref = Rand_ref+1;
PAA_c = 530000 + 20000*randn;       % PAA feed concentration [mg/L]

rng(Seed_ref + Batch_no + Rand_ref);
N_conc_paa = 2*75000 + 2000*randn;  % N in PAA feed [mg/L]

%% --- Process Disturbances ---
rng(Random_seed_ref + Batch_no);
b1 = 1 - 0.995;
a1 = [1 -0.995];

v = randn(T/h+1, 1);
distMuP = filter(b1,a1,0.03*v);
Xinterp.distMuP = createChannel('Penicillin specific growth rate disturbance','g/Lh','h',Batch_time,distMuP);

v = randn(T/h+1, 1);
distMuX = filter(b1,a1,0.25*v);
Xinterp.distMuX = createChannel('Biomass specific growth rate disturbance','hr^{-1}','h',Batch_time,distMuX);

v = randn(T/h+1, 1);
distcs = filter(b1,a1,5*300*v);
Xinterp.distcs = createChannel('Substrate concentration disturbance','g L^{-1}','h',Batch_time,distcs);

v = randn(T/h+1, 1);
distcoil = filter(b1,a1,300*v);
Xinterp.distcoil = createChannel('Oil concentration disturbance','g L^{-1}','h',Batch_time,distcoil);

v = randn(T/h+1, 1);
distabc = filter(b1,a1,0.2*v);
Xinterp.distabc = createChannel('Acid/Base concentration disturbance','g L^{-1}','h',Batch_time,distabc);

v = randn(T/h+1, 1);
distPAA = filter(b1,a1,300000*v);
Xinterp.distPAA = createChannel('PAA concentration disturbance','g L^{-1}','h',Batch_time,distPAA);

v = randn(T/h+1, 1);
distTcin = filter(b1,a1,100*v);
Xinterp.distTcin = createChannel('Coolant inlet temperature disturbance','K','h',Batch_time,distTcin);

v = randn(T/h+1, 1);
distO_2in = filter(b1,a1,0.02*v);
Xinterp.distO_2in = createChannel('Oxygen inlet concentration','%','h',Batch_time,distO_2in);

%% --- Load Parameter List ---
par = Parameter_list(x0, alpha_kla, N_conc_paa, PAA_c);

%% --- RUN SIMULATION ---
disp('Running IndPenSim... please wait (1-2 mins)');
Xref = indpensim(@fctrl_indpensim, Xinterp, x0, h, T, 2, par, Ctrl_flags);
disp('Simulation complete!');

%% --- Extract Variables ---
time       = Xref.P.t;       % time vector (hours)
penicillin = Xref.P.y;       % penicillin concentration [g/L]
biomass    = Xref.X.y;       % biomass [g/L]
substrate  = Xref.S.y;       % substrate [g/L]
dissolvedO2= Xref.DO2.y;     % dissolved oxygen [mg/L]
pH         = Xref.pH.y;      % pH
temperature= Xref.T.y;       % temperature [K]

%% --- PLOTS ---
figure('Name','IndPenSim Batch 1 - Baseline', 'NumberTitle','off', ...
       'Position',[100 100 1200 700]);

subplot(2,3,1)
plot(time, biomass, 'b', 'LineWidth', 1.5)
xlabel('Time (h)'); ylabel('g/L'); title('Biomass (X)'); grid on;

subplot(2,3,2)
plot(time, substrate, 'r', 'LineWidth', 1.5)
xlabel('Time (h)'); ylabel('g/L'); title('Substrate (S)'); grid on;

subplot(2,3,3)
plot(time, penicillin, 'g', 'LineWidth', 1.5)
xlabel('Time (h)'); ylabel('g/L'); title('Penicillin (P)'); grid on;

subplot(2,3,4)
plot(time, dissolvedO2, 'm', 'LineWidth', 1.5)
xlabel('Time (h)'); ylabel('mg/L'); title('Dissolved Oxygen (DO2)'); grid on;

subplot(2,3,5)
plot(time, pH, 'k', 'LineWidth', 1.5)
xlabel('Time (h)'); ylabel('pH'); title('pH'); grid on;

subplot(2,3,6)
plot(time, temperature - 273.15, 'c', 'LineWidth', 1.5)
xlabel('Time (h)'); ylabel('°C'); title('Temperature'); grid on;

sgtitle('IndPenSim — Batch 1 Baseline', 'FontSize', 14, 'FontWeight', 'bold')
saveas(gcf, '../results/batch1_overview.png')
disp('Plot saved to results folder')

%% --- Batch Summary ---
fprintf('\n========== BATCH SUMMARY ==========\n')
fprintf('Final Penicillin        : %.4f g/L\n', penicillin(end))
fprintf('Peak Biomass            : %.4f g/L\n', max(biomass))
fprintf('Min Dissolved Oxygen    : %.4f mg/L\n', min(dissolvedO2))
fprintf('pH Range                : %.2f - %.2f\n', min(pH), max(pH))
fprintf('Temp Range              : %.2f - %.2f C\n', ...
        min(temperature)-273.15, max(temperature)-273.15)
% Penicillin yield stats (calculated manually)
Pen_harvested = sum(Xref.Fremoved.y .* Xref.P.y) * h;
Pen_end       = Xref.V.y(end) * Xref.P.y(end);
Pen_total     = Pen_end - Pen_harvested;

fprintf('\n--- Yield Stats ---\n')
fprintf('Penicillin harvested during batch : %d Kg\n', round(Pen_harvested/1000))
fprintf('Final penicillin at harvest       : %d Kg\n', round(Pen_end/1000))
fprintf('Total penicillin yield            : %d Kg\n', round(Pen_total/1000))
fprintf('====================================\n')