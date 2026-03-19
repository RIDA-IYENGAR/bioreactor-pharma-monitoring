clc; clear; close all;

%% --- Settings ---
h          = 0.2;   % sampling time (hours)
T          = 230;   % batch duration (hours)
N_batches  = 20;    % number of batches to run

% Pre-allocate storage (T/h points, not T/h+1)
n_steps = T/h;
all_penicillin  = zeros(n_steps, N_batches);
all_biomass     = zeros(n_steps, N_batches);
all_substrate   = zeros(n_steps, N_batches);
all_DO2         = zeros(n_steps, N_batches);
all_pH          = zeros(n_steps, N_batches);
all_temp        = zeros(n_steps, N_batches);

% Storage for final yield stats
final_pen    = zeros(N_batches, 1);
peak_biomass = zeros(N_batches, 1);
pen_yield_kg = zeros(N_batches, 1);

fprintf('Running %d batches...\n', N_batches);

%% --- Run all batches ---
for Batch_no = 1:N_batches

    fprintf('  Batch %d / %d\n', Batch_no, N_batches);

    Batch_time = 0:h:T;

    % Control flags
    Ctrl_flags.SBC        = 0;
    Ctrl_flags.PRBS       = 0;
    Ctrl_flags.Fixed_Batch_length = 0;
    Ctrl_flags.IC         = 0;
    Ctrl_flags.Inhib      = 2;
    Ctrl_flags.Dis        = 1;  % disturbances ON - creates batch variation
    Ctrl_flags.Faults     = 0;
    Ctrl_flags.Vis        = 0;
    Ctrl_flags.Raman_spec = 0;
    Ctrl_flags.Batch_Num  = Batch_no;
    Ctrl_flags.Off_line_m     = 12;
    Ctrl_flags.Off_line_delay =  4;
    Ctrl_flags.plots      = 0;
    Ctrl_flags.T_sp       = 298;
    Ctrl_flags.pH_sp      = 6.5;

    % Each batch gets a different random seed = batch-to-batch variation
    Random_seed_ref = Batch_no * 10;
    Seed_ref = 31 + Random_seed_ref;
    Rand_ref = 1;

    % Initial conditions
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
    x0.a3             = 0;
    x0.a4             = 0;
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
    b1 = 1-0.995; a1 = [1 -0.995];

    v=randn(T/h+1,1); distMuP=filter(b1,a1,0.03*v);
    Xinterp.distMuP=createChannel('Pen growth disturbance','g/Lh','h',Batch_time,distMuP);
    v=randn(T/h+1,1); distMuX=filter(b1,a1,0.25*v);
    Xinterp.distMuX=createChannel('Biomass growth disturbance','hr^{-1}','h',Batch_time,distMuX);
    v=randn(T/h+1,1); distcs=filter(b1,a1,5*300*v);
    Xinterp.distcs=createChannel('Substrate disturbance','g L^{-1}','h',Batch_time,distcs);
    v=randn(T/h+1,1); distcoil=filter(b1,a1,300*v);
    Xinterp.distcoil=createChannel('Oil disturbance','g L^{-1}','h',Batch_time,distcoil);
    v=randn(T/h+1,1); distabc=filter(b1,a1,0.2*v);
    Xinterp.distabc=createChannel('Acid/Base disturbance','g L^{-1}','h',Batch_time,distabc);
    v=randn(T/h+1,1); distPAA=filter(b1,a1,300000*v);
    Xinterp.distPAA=createChannel('PAA disturbance','g L^{-1}','h',Batch_time,distPAA);
    v=randn(T/h+1,1); distTcin=filter(b1,a1,100*v);
    Xinterp.distTcin=createChannel('Coolant disturbance','K','h',Batch_time,distTcin);
    v=randn(T/h+1,1); distO_2in=filter(b1,a1,0.02*v);
    Xinterp.distO_2in=createChannel('O2 inlet disturbance','%','h',Batch_time,distO_2in);

    % Parameter list and run
    par  = Parameter_list(x0, alpha_kla, N_conc_paa, PAA_c);
    Xref = indpensim(@fctrl_indpensim, Xinterp, x0, h, T, 2, par, Ctrl_flags);

    % Store results
    all_penicillin(:,Batch_no) = Xref.P.y;
    all_biomass(:,Batch_no)    = Xref.X.y;
    all_substrate(:,Batch_no)  = Xref.S.y;
    all_DO2(:,Batch_no)        = Xref.DO2.y;
    all_pH(:,Batch_no)         = Xref.pH.y;
    all_temp(:,Batch_no)       = Xref.T.y - 273.15;

    % Yield stats
    final_pen(Batch_no)    = Xref.P.y(end);
    peak_biomass(Batch_no) = max(Xref.X.y);
    pen_yield_kg(Batch_no) = (Xref.V.y(end) * Xref.P.y(end)) / 1000;

end

time = Xref.P.t;
fprintf('All batches complete!\n');

%% --- PLOT 1: Overlay all 20 batches ---
figure('Name','20 Batch Overlay','NumberTitle','off','Position',[50 50 1300 750]);
colors = lines(N_batches); % 20 distinct colours

vars    = {all_penicillin, all_biomass, all_substrate, all_DO2, all_pH, all_temp};
titles  = {'Penicillin (P)', 'Biomass (X)', 'Substrate (S)', ...
           'Dissolved O2', 'pH', 'Temperature (°C)'};
ylabels = {'g/L','g/L','g/L','mg/L','pH','°C'};

for v = 1:6
    subplot(2,3,v)
    hold on
    for b = 1:N_batches
        plot(time, vars{v}(:,b), 'Color', [colors(b,:) 0.5], 'LineWidth', 1)
    end
    xlabel('Time (h)'); ylabel(ylabels{v});
    title(titles{v}); grid on; box on;
end
sgtitle('IndPenSim — 20 Batch Overlay (Natural Variation)', ...
        'FontSize',14,'FontWeight','bold')
saveas(gcf, '../results/20batch_overlay.png')

%% --- PLOT 2: Final penicillin yield distribution ---
figure('Name','Penicillin Yield Distribution','NumberTitle','off', ...
       'Position',[100 100 600 400]);
bar(1:N_batches, final_pen, 'FaceColor',[0.2 0.6 0.4])
hold on
yline(mean(final_pen), 'r--', 'LineWidth', 2, 'Label', ...
      sprintf('Mean = %.2f g/L', mean(final_pen)))
xlabel('Batch Number'); ylabel('Final Penicillin (g/L)')
title('Batch-to-Batch Variation in Penicillin Yield')
grid on
saveas(gcf, '../results/penicillin_yield_distribution.png')

%% --- Summary Statistics ---
fprintf('\n====== 20 BATCH STATISTICS ======\n')
fprintf('Penicillin  — Mean: %.2f  SD: %.2f  Min: %.2f  Max: %.2f g/L\n', ...
        mean(final_pen), std(final_pen), min(final_pen), max(final_pen))
fprintf('Peak Biomass— Mean: %.2f  SD: %.2f g/L\n', ...
        mean(peak_biomass), std(peak_biomass))
fprintf('Yield (Kg)  — Mean: %.0f  SD: %.0f Kg\n', ...
        mean(pen_yield_kg), std(pen_yield_kg))
fprintf('=================================\n')

%% --- Save data for Phase 2 ---
save('../results/batch_data.mat', ...
     'all_penicillin','all_biomass','all_substrate', ...
     'all_DO2','all_pH','all_temp','time', ...
     'final_pen','peak_biomass','pen_yield_kg')
fprintf('Batch data saved to results/batch_data.mat\n')