 
%% =========================================================
% phase5_substrate_PID_controller.m
% Phase 5: PID Controller for Substrate Feed Rate
%
% OBJECTIVE:
% Design a PID controller that maintains substrate
% concentration stable throughout the batch — identified
% as the TOP CPP in Phase 2 analysis.
% Then integrate with Phase 4 pH controller for full
% dual-loop control.
%
% CONTRIBUTOR: [Collaborator 2 name]
% =========================================================
clc; clear; close all;

%% =========================================================
%  SECTION 1: Setup — same as Phase 4
% =========================================================

fprintf('Setting up simulation...\n');

h = 0.2; T = 230; Batch_no = 1;
Batch_time = 0:h:T;

Ctrl_flags.SBC=0; Ctrl_flags.PRBS=0;
Ctrl_flags.Fixed_Batch_length=0; Ctrl_flags.IC=0;
Ctrl_flags.Inhib=2; Ctrl_flags.Dis=1; Ctrl_flags.Faults=0;
Ctrl_flags.Vis=0; Ctrl_flags.Raman_spec=0;
Ctrl_flags.Batch_Num=Batch_no; Ctrl_flags.Off_line_m=12;
Ctrl_flags.Off_line_delay=4; Ctrl_flags.plots=0;
Ctrl_flags.T_sp=298; Ctrl_flags.pH_sp=6.5;

Random_seed_ref=42; Seed_ref=31+Random_seed_ref; Rand_ref=1;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1; intial_conds=0.5+0.05*randn;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1; x0.mux=0.41+0.025*randn;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1; x0.mup=0.041+0.0025*randn;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1; x0.S=1+0.1*randn;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1; x0.DO2=15+0.5*randn;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1; x0.X=intial_conds+0.1*randn;
x0.P=0;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1; x0.V=5.800e+04+500*randn;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1; x0.Wt=6.2e+04+500*randn;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1; x0.CO2outgas=0.038+0.001*randn;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1; x0.O2=0.20+0.05*randn;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1; x0.pH=6.5+0.1*randn;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1; x0.T=297+0.5*randn;
x0.a0=intial_conds*(1/3); x0.a1=intial_conds*(2/3);
x0.a3=0; x0.a4=0; x0.Culture_age=0;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1; x0.PAA=1400+50*randn;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1; x0.NH3=1700+50*randn;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1; alpha_kla=85+10*randn;
rng(Seed_ref+Batch_no+Rand_ref); Rand_ref=Rand_ref+1; PAA_c=530000+20000*randn;
rng(Seed_ref+Batch_no+Rand_ref); N_conc_paa=2*75000+2000*randn;

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

% Run baseline first
fprintf('Running baseline batch...\n');
Xref_baseline = indpensim(@fctrl_indpensim, Xinterp, x0, h, T, 2, par, Ctrl_flags);
fprintf('Baseline complete.\n');

%% =========================================================
%  SECTION 2: Substrate PID Controller
%
%  Substrate (glucose) is the TOP CPP from Phase 2.
%  Target: keep substrate near 0.5 g/L throughout batch
%  Control output: glucose feed rate (Fs)
%
%  TUNING GUIDE:
%  Start with Kp=0.1, Ki=0.005, Kd=0.0005
%  If substrate oscillates → reduce Kp
%  If substrate drifts → increase Ki
% =========================================================

% PID gains for substrate control - TUNE THESE
Kp_s = 0.1;
Ki_s = 0.005;
Kd_s = 0.0005;

S_setpoint = 0.5;  % target substrate concentration g/L

% PID state
pid_s.integral  = 0;
pid_s.prev_error = 0;
pid_s.Kp = Kp_s;
pid_s.Ki = Ki_s;
pid_s.Kd = Kd_s;
pid_s.setpoint = S_setpoint;

    function [u, Xout] = substrate_PID_controller(X, Xd, k, h_step, T_end, Ctrl_flags_in)
        % Read current substrate
        if k > 1 && isfield(X,'S') && length(X.S.y) >= k
            current_S = X.S.y(k);
        else
            current_S = pid_s.setpoint;
        end

        % PID calculation
        error = pid_s.setpoint - current_S;
        pid_s.integral   = pid_s.integral + error * h_step;
        derivative       = (error - pid_s.prev_error) / h_step;
        pid_s.prev_error = error;

        u_pid = pid_s.Kp * error + ...
                pid_s.Ki * pid_s.integral + ...
                pid_s.Kd * derivative;

        % Get default recipe inputs
        u = fctrl_indpensim(X, Xd, k, h_step, T_end);

        % Apply PID to glucose feed rate
        % Clamp between 0 and 0.1 L/h (physical limits)
        u.Fs = max(0, min(0.1, u.Fs + u_pid));
        Xout = X;
    end

fprintf('Running substrate PID controlled batch...\n');
Xref_sub = indpensim(@substrate_PID_controller, Xinterp, x0, h, T, 2, par, Ctrl_flags);
fprintf('Substrate PID batch complete.\n');

%% =========================================================
%  SECTION 3: Dual Loop — pH + Substrate control together
%  This is the full integrated controller
% =========================================================

% Reset PID states for dual loop run
pid_s.integral=0; pid_s.prev_error=0;
pid_pH.integral=0; pid_pH.prev_error=0;
pid_pH.Kp=0.5; pid_pH.Ki=0.01; pid_pH.Kd=0.001;
pid_pH.setpoint=6.5;

    function u = dual_PID_controller(X, Xd, k, h_step, T_end, Ctrl_flags_in)
        u = fctrl_indpensim(X, Xd, k, h_step, T_end, Ctrl_flags_in);

        % --- Substrate PID ---
        if k>1 && isfield(X,'S') && length(X.S.y)>=k
            current_S = X.S.y(k);
        else
            current_S = pid_s.setpoint;
        end
        err_s = pid_s.setpoint - current_S;
        pid_s.integral = pid_s.integral + err_s*h_step;
        deriv_s = (err_s - pid_s.prev_error)/h_step;
        pid_s.prev_error = err_s;
        u_s = pid_s.Kp*err_s + pid_s.Ki*pid_s.integral + pid_s.Kd*deriv_s;
        u.Fs = max(0, min(0.1, u.Fs + u_s));

        % --- pH PID ---
        if k>1 && isfield(X,'pH') && length(X.pH.y)>=k
            current_pH = X.pH.y(k);
        else
            current_pH = pid_pH.setpoint;
        end
        err_pH = pid_pH.setpoint - current_pH;
        pid_pH.integral = pid_pH.integral + err_pH*h_step;
        deriv_pH = (err_pH - pid_pH.prev_error)/h_step;
        pid_pH.prev_error = err_pH;
        u_pH = pid_pH.Kp*err_pH + pid_pH.Ki*pid_pH.integral + pid_pH.Kd*deriv_pH;
        u.Fb = max(0, min(0.05, u.Fb + u_pH));
    end

fprintf('Running dual-loop PID batch (pH + Substrate)...\n');
Xref_dual = indpensim(@dual_PID_controller, Xinterp, x0, h, T, 2, par, Ctrl_flags);
fprintf('Dual loop batch complete.\n');

%% =========================================================
%  SECTION 4: Results and Plots
% =========================================================

time = Xref_baseline.S.t;

% --- Plot 1: Substrate control ---
figure('Name','Substrate Control','NumberTitle','off',...
       'Position',[100 100 1000 400]);
subplot(1,3,1)
plot(time, Xref_baseline.S.y,'b','LineWidth',1.5)
yline(S_setpoint,'r--','LineWidth',2,'Label','Setpoint')
xlabel('Time (h)'); ylabel('Substrate (g/L)')
title('Baseline'); grid on

subplot(1,3,2)
plot(time, Xref_sub.S.y,'g','LineWidth',1.5)
yline(S_setpoint,'r--','LineWidth',2,'Label','Setpoint')
xlabel('Time (h)'); ylabel('Substrate (g/L)')
title('Substrate PID Only'); grid on

subplot(1,3,3)
plot(time, Xref_dual.S.y,'m','LineWidth',1.5)
yline(S_setpoint,'r--','LineWidth',2,'Label','Setpoint')
xlabel('Time (h)'); ylabel('Substrate (g/L)')
title('Dual Loop (pH + Substrate)'); grid on

sgtitle('Substrate Control Comparison','FontSize',13,'FontWeight','bold')
saveas(gcf,'substrate_control_comparison.png')

% --- Plot 2: Penicillin yield all three ---
figure('Name','Penicillin Yield Comparison','NumberTitle','off',...
       'Position',[100 100 700 450]);
plot(time, Xref_baseline.P.y,'b-','LineWidth',2,'DisplayName','Baseline')
hold on
plot(time, Xref_sub.P.y,'g-','LineWidth',2,'DisplayName','Substrate PID')
plot(time, Xref_dual.P.y,'m-','LineWidth',2,'DisplayName','Dual Loop')
xlabel('Time (h)'); ylabel('Penicillin (g/L)')
title('Penicillin Yield: Baseline vs Controlled Batches')
legend('Location','northwest'); grid on
saveas(gcf,'penicillin_yield_all_controllers.png')

% --- Performance Summary ---
S_baseline_sd  = std(Xref_baseline.S.y);
S_sub_sd       = std(Xref_sub.S.y);
S_dual_sd      = std(Xref_dual.S.y);

fprintf('\n====== SUBSTRATE CONTROLLER PERFORMANCE ======\n')
fprintf('Substrate SD - Baseline    : %.4f g/L\n', S_baseline_sd)
fprintf('Substrate SD - Sub PID     : %.4f g/L\n', S_sub_sd)
fprintf('Substrate SD - Dual Loop   : %.4f g/L\n', S_dual_sd)
fprintf('Improvement (dual vs base) : %.1f%%\n', ...
        (1-S_dual_sd/S_baseline_sd)*100)
fprintf('\nFinal Penicillin:\n')
fprintf('  Baseline   : %.4f g/L\n', Xref_baseline.P.y(end))
fprintf('  Sub PID    : %.4f g/L\n', Xref_sub.P.y(end))
fprintf('  Dual Loop  : %.4f g/L\n', Xref_dual.P.y(end))
fprintf('===============================================\n')

% Save
save('phase5_results.mat','Xref_baseline','Xref_sub','Xref_dual',...
     'Kp_s','Ki_s','Kd_s','S_baseline_sd','S_sub_sd','S_dual_sd')
fprintf('Phase 5 results saved.\n') 
