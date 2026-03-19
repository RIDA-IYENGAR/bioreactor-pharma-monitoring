%% =========================================================
% phase2_CPP_analysis.m
% Phase 2: PCA and PLS to identify CPPs and CQAs
% =========================================================
clc; clear; close all;

%% --- Load batch data from Phase 1 ---
load('../results/batch_data.mat')
fprintf('Loaded batch data: %d batches\n', size(all_penicillin,2));

%% =========================================================
%  SECTION 1: Build feature matrix
%  For each batch, extract summary statistics as features
%  These represent what you'd know about a batch
% =========================================================

N_batches = size(all_penicillin, 2);

% For each batch, compute mean and SD of each variable
% This gives us a 20 x 12 feature matrix (6 vars x 2 stats)
features = zeros(N_batches, 12);
feat_names = {'DO2_mean','DO2_sd', ...
              'pH_mean','pH_sd', ...
              'Temp_mean','Temp_sd', ...
              'Biomass_mean','Biomass_sd', ...
              'Substrate_mean','Substrate_sd', ...
              'Pen_mean','Pen_sd'};

for b = 1:N_batches
    features(b,1)  = mean(all_DO2(:,b));
    features(b,2)  = std(all_DO2(:,b));
    features(b,3)  = mean(all_pH(:,b));
    features(b,4)  = std(all_pH(:,b));
    features(b,5)  = mean(all_temp(:,b));
    features(b,6)  = std(all_temp(:,b));
    features(b,7)  = mean(all_biomass(:,b));
    features(b,8)  = std(all_biomass(:,b));
    features(b,9)  = mean(all_substrate(:,b));
    features(b,10) = std(all_substrate(:,b));
    features(b,11) = mean(all_penicillin(:,b));
    features(b,12) = std(all_penicillin(:,b));
end

% Target variable: final penicillin concentration (our CQA)
Y = final_pen;

fprintf('Feature matrix built: %d batches x %d features\n', ...
        size(features,1), size(features,2));

%% =========================================================
%  SECTION 2: PCA
%  Find which combinations of variables explain most variance
% =========================================================

% Standardise features (zero mean, unit variance) - required for PCA
X_std = zscore(features);

% Run PCA
[coeff, score, latent, ~, explained] = pca(X_std);

% --- Plot 1: Scree plot (how much variance each PC explains) ---
figure('Name','PCA - Scree Plot','NumberTitle','off',...
       'Position',[100 100 600 400]);
bar(explained, 'FaceColor',[0.2 0.5 0.8])
hold on
plot(cumsum(explained), 'ro-', 'LineWidth', 2, 'MarkerSize', 8)
xlabel('Principal Component')
ylabel('Variance Explained (%)')
title('PCA Scree Plot — Variance Explained per Component')
legend('Individual','Cumulative','Location','east')
yline(80, 'k--', '80% threshold')
grid on
saveas(gcf,'../results/PCA_scree.png')

fprintf('\n--- PCA Results ---\n')
fprintf('PC1 explains: %.1f%%\n', explained(1))
fprintf('PC2 explains: %.1f%%\n', explained(2))
fprintf('PC1+PC2 total: %.1f%%\n', explained(1)+explained(2))

% --- Plot 2: PCA scores plot (each dot = one batch) ---
figure('Name','PCA - Scores Plot','NumberTitle','off',...
       'Position',[100 100 650 500]);

% Colour batches by final penicillin yield
scatter(score(:,1), score(:,2), 120, Y, 'filled')
colorbar; colormap(jet)
xlabel(sprintf('PC1 (%.1f%%)', explained(1)))
ylabel(sprintf('PC2 (%.1f%%)', explained(2)))
title('PCA Scores — Batches Coloured by Penicillin Yield')
grid on; box on

% Label each batch number
for b = 1:N_batches
    text(score(b,1)+0.05, score(b,2)+0.05, num2str(b), 'FontSize', 8)
end
saveas(gcf,'../results/PCA_scores.png')

% --- Plot 3: PCA loadings (which variables drive each PC) ---
figure('Name','PCA - Loadings','NumberTitle','off',...
       'Position',[100 100 800 400]);
bar(coeff(:,1:2))
set(gca,'XTickLabel', feat_names,'XTickLabelRotation',45)
legend('PC1','PC2')
title('PCA Loadings — Variable Contributions to PC1 and PC2')
ylabel('Loading value')
grid on
saveas(gcf,'../results/PCA_loadings.png')

%% =========================================================
%  SECTION 3: PLS Regression
%  Directly link process variables to penicillin yield
%  This identifies your CPPs
% =========================================================

% Use only process variables as predictors (not penicillin itself)
% Columns 1-10: DO2, pH, Temp, Biomass, Substrate (mean + SD each)
X_process = features(:, 1:10);
X_process_std = zscore(X_process);
Y_std = zscore(Y);

process_names = feat_names(1:10);

% Run PLS with cross-validation to find optimal number of components
[~,~,~,~,BETA,PCTVAR,~,stats] = plsregress(X_process_std, Y_std, 5);

% --- Plot 4: PLS variance explained ---
figure('Name','PLS - Variance Explained','NumberTitle','off',...
       'Position',[100 100 600 400]);
plot(1:5, cumsum(PCTVAR(2,:))*100, 'bo-', 'LineWidth', 2, 'MarkerSize', 8)
xlabel('Number of PLS Components')
ylabel('Cumulative Variance Explained in Y (%)')
title('PLS — Variance in Penicillin Yield Explained')
yline(80,'k--','80% threshold')
grid on
saveas(gcf,'../results/PLS_variance.png')

% --- Plot 5: PLS coefficients = CPP importance ranking ---
figure('Name','PLS - CPP Importance','NumberTitle','off',...
       'Position',[100 100 800 450]);
beta_vals = BETA(2:end); % exclude intercept
[sorted_beta, sort_idx] = sort(abs(beta_vals), 'descend');

barh(sorted_beta, 'FaceColor',[0.8 0.3 0.3])
set(gca,'YTickLabel', process_names(sort_idx))
xlabel('|PLS Coefficient| — Importance for Penicillin Yield')
title('CPP Ranking: Process Variables vs Penicillin Yield')
grid on
saveas(gcf,'../results/PLS_CPP_ranking.png')

% --- Plot 6: PLS predicted vs actual yield ---
% Manually compute predictions using PLS coefficients
[~,~,~,~,BETA3] = plsregress(X_process_std, Y_std, 3);
Ypred_std = [ones(size(X_process_std,1),1) X_process_std] * BETA3;
Ypred = Ypred_std * std(Y) + mean(Y);

figure('Name','PLS - Predicted vs Actual','NumberTitle','off',...
       'Position',[100 100 550 500]);
scatter(Y, Ypred, 100, 'b', 'filled')
hold on
plot([min(Y) max(Y)], [min(Y) max(Y)], 'r--', 'LineWidth', 2)
xlabel('Actual Final Penicillin (g/L)')
ylabel('PLS Predicted Penicillin (g/L)')
title('PLS Model: Predicted vs Actual Penicillin Yield')
grid on; box on

% R-squared
SS_res = sum((Y - Ypred).^2);
SS_tot = sum((Y - mean(Y)).^2);
R2 = 1 - SS_res/SS_tot;
text(min(Y)+1, max(Y)-2, sprintf('R² = %.3f', R2), 'FontSize', 12, ...
     'FontWeight','bold','Color','r')
saveas(gcf,'../results/PLS_predicted_vs_actual.png')

%% --- Print CPP Summary ---
fprintf('\n====== CPP IDENTIFICATION RESULTS ======\n')
fprintf('Ranked by PLS coefficient magnitude:\n\n')
for i = 1:length(sort_idx)
    fprintf('  %d. %-20s  |beta| = %.4f\n', ...
            i, process_names{sort_idx(i)}, sorted_beta(i))
end
fprintf('\nPLS Model R² = %.3f\n', R2)
fprintf('=========================================\n')

%% --- Save results ---
save('../results/phase2_results.mat', ...
     'features','feat_names','coeff','score','explained', ...
     'BETA','R2','Y','Ypred','process_names')
fprintf('Phase 2 results saved.\n')