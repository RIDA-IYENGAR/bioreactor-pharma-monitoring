%% =========================================================
% phase4_pH_PID_controller.m
% Phase 4: PID Controller for pH in IndPenSim
%
% OBJECTIVE:
% Design and tune a PID controller that maintains pH at
% setpoint (6.5) throughout the batch, reducing the pH
% variation identified as a top CPP in Phase 2.
%
% CONTRIBUTOR: [Collaborator 1 name]
% =========================================================
clc; clear; close all;

%% =========================================================
%  SECTION 1: Run uncontrolled batch (baseline for comparison)
%  This shows what pH looks like WITHOUT your controller
% =========================================================

fprintf('Running baseline batch (no pH control)...\n');

h = 0.2;
T = 230;
Batch_no = 1;
Batch_time = 0:h:T;

% Control flags - standard run
Ctrl_flags.SBC        = 0;
Ctrl_flags.PRBS       = 0;
Ctrl_flags.Fixed_Batch_length = 0;
Ctrl_flags.IC         = 0;
Ctrl_flags.Inhib      = 2;
Ctrl_flags.Dis        = 1;
Ctrl_flags.Faults     = 0;
Ctrl_flags.Vis        = 0;
Ctrl_flags.Raman_spec = 0;
Ctrl_flags.Batch_Num  = Batch_no;
Ctrl_flags.Off_line_m     = 12;
Ctrl_flags.Off_line_delay =  4;
Ctrl_flags.plots      = 0;
Ctrl_flags.T_sp       = 298;
Ctrl_flags.pH_sp      = 6.5;  % setpoint - controller must track this

% Initial conditions
Random_seed_ref = 42;
Seed_ref = 31 + Random_seed_ref;
Rand_ref = 1;

rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1;
intial_conds      = 0.5+0.05*randn;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1;
x0.mux            = 0.41+0.025*randn;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1;
x0.mup            = 0.041+0.0025*randn;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1;
x0.S              = 1+0.1*randn;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1;
x0.DO2            = 15+0.5*randn;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1;
x0.X              = intial_conds+0.1*randn;
x0.P              = 0;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1;
x0.V              = 5.800e+04+500*randn;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1;
x0.Wt             = 6.2e+04+500*randn;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1;
x0.CO2outgas      = 0.038+0.001*randn;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1;
x0.O2             = 0.20+0.05*randn;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1;
x0.pH             = 6.5+0.1*randn;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1;
x0.T              = 297+0.5*randn;
x0.a0             = intial_conds*(1/3);
x0.a1             = intial_conds*(2/3);
x0.a3             = 0; x0.a4 = 0;
x0.Culture_age    = 0;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1;
x0.PAA            = 1400+50*randn;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1;
x0.NH3            = 1700+50*randn;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1;
alpha_kla         = 85+10*randn;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1;
PAA_c             = 530000+20000*randn;
rng(Seed_ref+Batch_no+Rand_ref);
N_conc_paa        = 2*75000+2000*randn;

% Disturbances
rng(Random_seed_ref+Batch_no);
b1=1-0.995; a1=[1 -0.995];
v=randn(T/h+1,1); Xinterp.distMuP=createChannel('distMuP','g/Lh','h',Batch_time,filter(b1,a1,0.03*v));
v=randn(T/h+1,1); Xinterp.distMuX=createChannel('distMuX','hr^{-1}','h',Batch_time,filter(b1,a1,0.25*v));
v=randn(T/h+1,1); Xinterp.distcs=createChannel('distcs','g L^{-1}','h',Batch_time,filter(b1,a1,5*300*v));
v=randn(T/h+1,1); Xinterp.distcoil=createChannel('distcoil','g L^{-1}','h',Batch_time,filter(b1,a1,300*v));
v=randn(T/h+1,1); Xinterp.distabc=createChannel('distabc','g L^{-1}','h',Batch_time,filter(b1,a1,0.2*v));
v=randn(T/h+1,1); Xinterp.distPAA=createChannel('distPAA','g L^{-1}','h',Batch_time,filter(b1,a1,300000*v));
v=randn(T/h+1,1); Xinterp.distTcin=createChannel('distTcin','K','h',Batch_time,filter(b1,a1,100*v));
v=randn(T/h+1,1); Xinterp.distO_2in=createChannel('distO2in','%','h',Batch_time,filter(b1,a1,0.02*v));

par = Parameter_list(x0, alpha_kla, N_conc_paa, PAA_c);
Xref_baseline = indpensim(@fctrl_indpensim, Xinterp, x0, h, T, 2, par, Ctrl_flags);
fprintf('Baseline complete.\n');

%% =========================================================
%  SECTION 2: PID Controller Implementation
%
%  HOW A PID CONTROLLER WORKS:
%  error = setpoint - measured_value
%  u(t) = Kp*error + Ki*integral(error) + Kd*derivative(error)
%
%  For pH control:
%  - Measured value = current pH
%  - Setpoint = 6.5
%  - Control output = base flow rate (Fbase) to correct pH
%
%  TUNING GUIDE:
%  Start with Kp=0.5, Ki=0.01, Kd=0.001
%  If pH oscillates → reduce Kp
%  If pH is slow to correct → increase Ki
%  If pH overshoots → increase Kd
% =========================================================

% PID Gains - TUNE THESE to get best performance
Kp = 0.5;    % Proportional gain
Ki = 0.01;   % Integral gain
Kd = 0.001;  % Derivative gain

pH_setpoint = 6.5;

% PID state variables
integral_error = 0;
prev_error     = 0;

% Storage for results
n_steps   = length(0:h:T);
pH_pid    = zeros(n_steps, 1);
Fbase_pid = zeros(n_steps, 1);
error_log = zeros(n_steps, 1);

fprintf('Running PID controlled batch...\n');

%% =========================================================
%  SECTION 3: Custom control function with PID
%  IndPenSim calls this function at every timestep
%  It reads current pH and adjusts base flow rate
% =========================================================

function u = pH_PID_controller(X, Xd, k, h_step, T_end, Ctrl_flags_in)

    % Persistent variables retain values between function calls
    persistent integral_error prev_error

    % Initialize persistent variables on first call
    if isempty(integral_error)
        integral_error = 0;
    end

    if isempty(prev_error)
        prev_error = 0;
    end

    % PID tuning parameters
    Kp = 0.5;
    Ki = 0.01;
    Kd = 0.001;

    % Setpoint
    pH_setpoint = 6.5;

    % Read current pH
    if k > 1 && isfield(X,'pH') && length(X.pH.y) >= k
        current_pH = X.pH.y(k);
    else
        current_pH = pH_setpoint;
    end

    % Calculate error
    error = pH_setpoint - current_pH;

    % Integral term
    integral_error = integral_error + error * h_step;

    % Derivative term
    derivative = (error - prev_error) / h_step;

    % Save error for next iteration
    prev_error = error;

    % PID equation
    u_pid = Kp*error + Ki*integral_error + Kd*derivative;

    % Default IndPenSim controls
    u = fctrl_indpensim(X, Xd, k, h_step, T_end, Ctrl_flags_in);
    % Apply PID correction to base flow
    u.Fb = max(0, min(0.05, u.Fb + u_pid));

end

% Run simulation with PID controller
Xref_pid = indpensim(@pH_PID_controller, Xinterp, x0, h, T, 2, par, Ctrl_flags);
fprintf('PID controlled batch complete.\n');

%% =========================================================
%  SECTION 4: Compare baseline vs PID controlled
% =========================================================

time = Xref_baseline.pH.t;

% --- Plot 1: pH comparison ---
figure('Name','pH Control - Baseline vs PID','NumberTitle','off',...
       'Position',[100 100 900 400]);

subplot(1,2,1)
plot(time, Xref_baseline.pH.y, 'b-', 'LineWidth', 1.5)
hold on
yline(pH_setpoint, 'r--', 'LineWidth', 2, 'Label', 'Setpoint 6.5')
xlabel('Time (h)'); ylabel('pH')
title('Baseline — No pH Control')
ylim([6.0 7.0]); grid on

subplot(1,2,2)
plot(time, Xref_pid.pH.y, 'g-', 'LineWidth', 1.5)
hold on
yline(pH_setpoint, 'r--', 'LineWidth', 2, 'Label', 'Setpoint 6.5')
xlabel('Time (h)'); ylabel('pH')
title(sprintf('PID Controlled — Kp=%.2f Ki=%.3f Kd=%.4f', Kp, Ki, Kd))
ylim([6.0 7.0]); grid on

sgtitle('pH Control: Baseline vs PID Controller','FontSize',13,'FontWeight','bold')
saveas(gcf,'pH_control_comparison.png')

% --- Plot 2: Penicillin comparison ---
figure('Name','Penicillin Yield - Baseline vs PID','NumberTitle','off',...
       'Position',[100 100 700 400]);
plot(time, Xref_baseline.P.y, 'b-', 'LineWidth', 2, 'DisplayName','Baseline')
hold on
plot(time, Xref_pid.P.y, 'g-', 'LineWidth', 2, 'DisplayName','PID Controlled')
xlabel('Time (h)'); ylabel('Penicillin (g/L)')
title('Penicillin Yield: Baseline vs pH PID Controlled')
legend('Location','northwest'); grid on
saveas(gcf,'pH_penicillin_comparison.png')

% --- Performance Summary ---
pH_baseline_sd = std(Xref_baseline.pH.y);
pH_pid_sd      = std(Xref_pid.pH.y);
pH_improvement = (1 - pH_pid_sd/pH_baseline_sd)*100;

fprintf('\n====== pH CONTROLLER PERFORMANCE ======\n')
fprintf('PID Gains: Kp=%.3f  Ki=%.4f  Kd=%.5f\n', Kp, Ki, Kd)
fprintf('pH SD - Baseline    : %.4f\n', pH_baseline_sd)
fprintf('pH SD - PID         : %.4f\n', pH_pid_sd)
fprintf('pH variation reduced: %.1f%%\n', pH_improvement)
fprintf('Final Pen - Baseline: %.4f g/L\n', Xref_baseline.P.y(end))
fprintf('Final Pen - PID     : %.4f g/L\n', Xref_pid.P.y(end))
fprintf('========================================\n')

% Save results
save('phase4_results.mat','Xref_baseline','Xref_pid',...
     'Kp','Ki','Kd','pH_baseline_sd','pH_pid_sd')
fprintf('Phase 4 results saved.\n')
