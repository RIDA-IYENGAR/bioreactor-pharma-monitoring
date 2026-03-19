%% =========================================================
% phase3_raman_batch.m
% Phase 3 Step 1: Run batch with Raman spectroscopy enabled
% =========================================================
clc; clear; close all;

%% --- Settings ---
h = 0.2;
T = 230;
N_batches = 10;  % run 10 batches with Raman for training data

fprintf('Running %d batches WITH Raman spectroscopy...\n', N_batches);

% Storage
raman_data   = {};   % cell array - each batch may have different length
pen_conc     = {};   % corresponding penicillin concentrations
biomass_conc = {};   % biomass concentrations

for Batch_no = 1:N_batches

    fprintf('  Batch %d / %d\n', Batch_no, N_batches);
    Batch_time = 0:h:T;

    % Control flags - Raman ON this time
    Ctrl_flags.SBC        = 0;
    Ctrl_flags.PRBS       = 0;
    Ctrl_flags.Fixed_Batch_length = 0;
    Ctrl_flags.IC         = 0;
    Ctrl_flags.Inhib      = 2;
    Ctrl_flags.Dis        = 1;
    Ctrl_flags.Faults     = 0;
    Ctrl_flags.Vis        = 0;
    Ctrl_flags.Raman_spec = 1;  % ← Raman spectroscopy ON
    Ctrl_flags.Batch_Num  = Batch_no;
    Ctrl_flags.Off_line_m     = 12;
    Ctrl_flags.Off_line_delay =  4;
    Ctrl_flags.plots      = 0;
    Ctrl_flags.T_sp       = 298;
    Ctrl_flags.pH_sp      = 6.5;

    % Random seed
    Random_seed_ref = Batch_no * 7;
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
    b1=1-0.995; a1=[1 -0.995];

    v=randn(T/h+1,1); Xinterp.distMuP=createChannel('distMuP','g/Lh','h',Batch_time,filter(b1,a1,0.03*v));
    v=randn(T/h+1,1); Xinterp.distMuX=createChannel('distMuX','hr^{-1}','h',Batch_time,filter(b1,a1,0.25*v));
    v=randn(T/h+1,1); Xinterp.distcs=createChannel('distcs','g L^{-1}','h',Batch_time,filter(b1,a1,5*300*v));
    v=randn(T/h+1,1); Xinterp.distcoil=createChannel('distcoil','g L^{-1}','h',Batch_time,filter(b1,a1,300*v));
    v=randn(T/h+1,1); Xinterp.distabc=createChannel('distabc','g L^{-1}','h',Batch_time,filter(b1,a1,0.2*v));
    v=randn(T/h+1,1); Xinterp.distPAA=createChannel('distPAA','g L^{-1}','h',Batch_time,filter(b1,a1,300000*v));
    v=randn(T/h+1,1); Xinterp.distTcin=createChannel('distTcin','K','h',Batch_time,filter(b1,a1,100*v));
    v=randn(T/h+1,1); Xinterp.distO_2in=createChannel('distO2in','%','h',Batch_time,filter(b1,a1,0.02*v));

    % Run
    par  = Parameter_list(x0, alpha_kla, N_conc_paa, PAA_c);
    Xref = indpensim(@fctrl_indpensim, Xinterp, x0, h, T, 2, par, Ctrl_flags);

    % Extract Raman spectra and concentrations
    % Raman data is stored in Xref.Raman_Spec
    if isfield(Xref, 'Raman_Spec')
        % Intensity is 2200 wavelengths x 1150 timepoints
        % Transpose to 1150 timepoints x 2200 wavelengths (standard format)
        raman_data{Batch_no}   = Xref.Raman_Spec.Intensity';
        pen_conc{Batch_no}     = Xref.P.y;
        biomass_conc{Batch_no} = Xref.X.y;
        fprintf('    Raman data: %d timepoints x %d wavelengths\n', ...
                size(Xref.Raman_Spec.Intensity,2), size(Xref.Raman_Spec.Intensity,1))
    else
        fprintf('    WARNING: No Raman field found for batch %d\n', Batch_no)
    end
end

fprintf('All Raman batches complete!\n')
wavelengths = Xref.Raman_Spec.Wavelength;
% Save
save('../results/raman_batch_data.mat', ...
     'raman_data','pen_conc','biomass_conc','wavelengths')
fprintf('Raman data saved.\n')