%% =========================================================
% phase3_soft_sensor.m
% Phase 3 Step 2: Build and validate Raman soft sensor
% Predicts penicillin concentration from Raman spectra
% =========================================================
clc; clear; close all;

%% --- Load Raman batch data ---
load('../results/raman_batch_data.mat')
fprintf('Loaded Raman data: %d batches\n', length(raman_data));

%% =========================================================
%  SECTION 1: Assemble training dataset
%  Stack all batches into one big matrix
%  X = spectra (rows = timepoints, cols = wavelengths)
%  Y = penicillin concentration at each timepoint
% =========================================================

X_all = [];
Y_pen = [];
Y_bio = [];
batch_labels = [];

for b = 1:length(raman_data)
    n_time = size(raman_data{b}, 1);

    % Make sure spectra and concentration vectors are same length
    n = min(n_time, length(pen_conc{b}));

    X_all = [X_all;  raman_data{b}(1:n, :)];
    Y_pen = [Y_pen;  pen_conc{b}(1:n)];
    Y_bio = [Y_bio;  biomass_conc{b}(1:n)];
    batch_labels = [batch_labels; b*ones(n,1)];
end

fprintf('Total dataset: %d timepoints x %d wavelengths\n', size(X_all));

%% =========================================================
%  SECTION 2: Preprocessing
%  Standardise spectra - critical for PLS to work correctly
% =========================================================

% Mean-centre and scale each wavelength
X_mean = mean(X_all);
X_std_val = std(X_all);
X_std_val(X_std_val == 0) = 1; % avoid division by zero

X_proc = (X_all - X_mean) ./ X_std_val;

% Standardise Y
Y_mean = mean(Y_pen);
Y_std_val = std(Y_pen);
Y_proc = (Y_pen - Y_mean) ./ Y_std_val;

fprintf('Preprocessing complete.\n')

%% =========================================================
%  SECTION 3: Train/Test Split
%  Batches 1-7 = training, Batches 8-10 = test
%  This is leave-batch-out validation - realistic for bioprocess
% =========================================================

train_idx = batch_labels <= 7;
test_idx  = batch_labels > 7;

X_train = X_proc(train_idx, :);
Y_train = Y_proc(train_idx);
X_test  = X_proc(test_idx, :);
Y_test  = Y_pen(test_idx);   % keep in original units for interpretability

fprintf('Training set: %d timepoints\n', sum(train_idx))
fprintf('Test set    : %d timepoints\n', sum(test_idx))

%% =========================================================
%  SECTION 4: PLS Soft Sensor
%  Find optimal number of components first
% =========================================================

fprintf('\nFinding optimal PLS components...\n')

max_comp = 15;
RMSECV   = zeros(max_comp, 1);

for n_comp = 1:max_comp
    [~,~,~,~,BETA_cv] = plsregress(X_train, Y_train, n_comp, ...
                                    'cv', 10);  % 10-fold cross validation
    Ypred_cv = [ones(size(X_train,1),1) X_train] * BETA_cv;
    % Denormalise
    Ypred_cv_real = Ypred_cv * Y_std_val + Y_mean;
    Y_train_real  = Y_train  * Y_std_val + Y_mean;
    RMSECV(n_comp) = sqrt(mean((Y_train_real - Ypred_cv_real).^2));
end

[~, opt_comp] = min(RMSECV);
fprintf('Optimal PLS components: %d\n', opt_comp)

% --- Plot 1: RMSECV vs components ---
figure('Name','PLS Components Selection','NumberTitle','off',...
       'Position',[100 100 600 400]);
plot(1:max_comp, RMSECV, 'bo-', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(opt_comp, RMSECV(opt_comp), 'r*', 'MarkerSize', 15, 'LineWidth', 2)
xlabel('Number of PLS Components')
ylabel('RMSECV (g/L)')
title('Cross-Validation: Optimal Number of PLS Components')
legend('RMSECV','Optimal','Location','northeast')
grid on
saveas(gcf,'../results/PLS_components_selection.png')

%% =========================================================
%  SECTION 5: Train final model and test on held-out batches
% =========================================================

fprintf('\nTraining final PLS model with %d components...\n', opt_comp)
[~,~,~,~,BETA_final] = plsregress(X_train, Y_train, opt_comp);

% Predict on test set
Ypred_test_std = [ones(size(X_test,1),1) X_test] * BETA_final;
Ypred_test = Ypred_test_std * Y_std_val + Y_mean;

% Predict on training set too (for comparison)
Ypred_train_std = [ones(size(X_train,1),1) X_train] * BETA_final;
Ypred_train = Ypred_train_std * Y_std_val + Y_mean;
Y_train_real = Y_train * Y_std_val + Y_mean;

% Performance metrics
RMSEP = sqrt(mean((Y_test - Ypred_test).^2));
RMSEC = sqrt(mean((Y_train_real - Ypred_train).^2));
SS_res = sum((Y_test - Ypred_test).^2);
SS_tot = sum((Y_test - mean(Y_test)).^2);
R2_test = 1 - SS_res/SS_tot;

fprintf('\n====== SOFT SENSOR PERFORMANCE ======\n')
fprintf('PLS Components used  : %d\n',   opt_comp)
fprintf('RMSEC (training)     : %.4f g/L\n', RMSEC)
fprintf('RMSEP (test)         : %.4f g/L\n', RMSEP)
fprintf('R² on test batches   : %.4f\n',  R2_test)
fprintf('======================================\n')

%% =========================================================
%  SECTION 6: Plots
% =========================================================

% --- Plot 2: Predicted vs Actual (test set) ---
figure('Name','Soft Sensor - Predicted vs Actual','NumberTitle','off',...
       'Position',[100 100 550 500]);
scatter(Y_test, Ypred_test, 30, 'b', 'filled', 'MarkerFaceAlpha', 0.4)
hold on
plot([min(Y_test) max(Y_test)],[min(Y_test) max(Y_test)],'r--','LineWidth',2)
xlabel('Actual Penicillin (g/L)')
ylabel('Predicted Penicillin (g/L)')
title('Soft Sensor: Predicted vs Actual (Test Batches 8-10)')
text(min(Y_test)+0.5, max(Y_test)-1, ...
     sprintf('R² = %.3f\nRMSEP = %.3f g/L', R2_test, RMSEP), ...
     'FontSize',11,'FontWeight','bold','Color','r')
grid on; box on
saveas(gcf,'../results/soft_sensor_predicted_vs_actual.png')

% --- Plot 3: Real-time prediction across one test batch ---
test_batch_times = find(batch_labels == 8);
time_vec = (0:length(test_batch_times)-1) * 0.2;

figure('Name','Soft Sensor - Real Time Prediction','NumberTitle','off',...
       'Position',[100 100 800 400]);
plot(time_vec, Y_test(1:length(test_batch_times)), ...
     'b-', 'LineWidth', 2, 'DisplayName', 'Actual')
hold on
plot(time_vec, Ypred_test(1:length(test_batch_times)), ...
     'r--', 'LineWidth', 2, 'DisplayName', 'Soft Sensor Predicted')
xlabel('Time (hours)')
ylabel('Penicillin Concentration (g/L)')
title('Soft Sensor: Real-Time Penicillin Prediction — Batch 8')
legend('Location','northwest')
grid on
saveas(gcf,'../results/soft_sensor_realtime.png')

% --- Plot 4: Regression coefficients (which wavelengths matter most) ---
beta_spectrum = BETA_final(2:end);  % exclude intercept
figure('Name','Soft Sensor - Important Wavelengths','NumberTitle','off',...
       'Position',[100 100 800 400]);
plot(wavelengths, beta_spectrum, 'k-', 'LineWidth', 1)
xlabel('Raman Wavelength (cm^{-1})')
ylabel('PLS Regression Coefficient')
title('Soft Sensor: Wavelengths Most Predictive of Penicillin')
grid on
saveas(gcf,'../results/soft_sensor_wavelengths.png')

%% --- Save model ---
save('../results/soft_sensor_model.mat', ...
     'BETA_final','X_mean','X_std_val','Y_mean','Y_std_val', ...
     'opt_comp','RMSEP','RMSEC','R2_test','wavelengths')
fprintf('Soft sensor model saved.\n')
