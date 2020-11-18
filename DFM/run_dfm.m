%% LP vs VAR: DFM SIMULATION STUDY
% Dake Li and Christian Wolf


%% HOUSEKEEPING

clear all
close all

addpath(genpath(fullfile('..', 'Auxiliary_Functions')))
addpath(genpath(fullfile('..', 'Estimation_Routines')))
addpath(genpath('Subroutines'))

rng(1, 'twister');
tic;

% Parallel computing object

num_partition = str2num(getenv('SLURM_ARRAY_TASK_COUNT'));
idx_partition = str2num(getenv('SLURM_ARRAY_TASK_ID'));
if isempty(num_partition)
    num_partition = 1;
    idx_partition = 1;
end
    
num_workers = str2num(getenv('SLURM_CPUS_PER_TASK'));
if ~isempty(num_workers)
    poolobj = parpool('local', num_workers);
else
    poolobj = parpool('local');
end
clear num_workers;


%% DECIDE WHICH EXPERIMENT TO RUN

% manually set up experiment

dgp_type = 'G'; % 'MP'; % Either 'G' or 'MP'
estimand_type = 'ObsShock'; % 'Recursive'; 'IV'; % Either 'ObsShock', 'Recursive', or 'IV'
lag_type = 4; % No. of lags to impose in estimation, or NaN (meaning AIC)

% overwrite the above experiment setup if specify another in bash

if ~isempty(getenv('MATLAB_DGP_TYPE'))
    dgp_type = getenv('MATLAB_DGP_TYPE');
    estimand_type = getenv('MATLAB_ESTIMAND_TYPE');
    lag_type = str2num(getenv('MATLAB_LAG_TYPE'));
end


%% SETTINGS

% Apply shared settings as well as settings specific to DGP and estimand type

run(fullfile('Settings', 'shared'));
run(fullfile('Settings', dgp_type));
run(fullfile('Settings', estimand_type));

% Storage folder for results

if isnan(lag_type)
    save_suff = '_aic';
else
    save_suff = num2str(lag_type);
end
save_folder = fullfile('Results', strcat('lag', save_suff));


%% DGP

%----------------------------------------------------------------
% Set up DGP
%----------------------------------------------------------------

% estimate DFM from dataset

DFM_estimate = DFM_est(DF_model.n_fac, DF_model.n_lags_fac, DF_model.n_lags_uar);

% store estimated DFM parameters

DF_model.Phi           = DFM_estimate.Phi;
DF_model.Sigma_eta     = DFM_estimate.Sigma_eta;

DF_model.Lambda        = DFM_estimate.Lambda;
DF_model.delta         = DFM_estimate.delta;
DF_model.sigma_v       = DFM_estimate.sigma_v;

DF_model.variable_name = DFM_estimate.bplabvec_long;

%----------------------------------------------------------------
% Calibrate IV strength and Shock Weight
%----------------------------------------------------------------

% extract factor shock series and external shock series

DFM_estimate.fac_shock                 = DFM_estimate.fac_res / chol(DFM_estimate.Sigma_eta);
external_shock_data                    = readtable(strcat('external_shock_series_', dgp_type, '.csv'));
DFM_estimate.external_shock            = external_shock_data{:,2};
DFM_estimate.external_shock_time_range = [external_shock_data{1,1}, external_shock_data{end,1}, 4];
clear external_shock_data;

% regress external shock series on factor shock series to calibrate

DFM_estimate.calibrate_out       = calibrateIV(DFM_estimate);
DF_model.calibrated_shock_weight = DFM_estimate.calibrate_out.weight;

if strcmp(estimand_type, 'IV')
    if DF_model.IV.IV_strength_calibrate==1
        DF_model.IV.alpha = DFM_estimate.calibrate_out.alpha;
        DF_model.IV.sigma_v = DFM_estimate.calibrate_out.sigma_v;
    else
        DF_model.IV.alpha = DF_model.IV.manual_alpha;
        DF_model.IV.sigma_v = DF_model.IV.manual_sigma_v;
    end
end

%----------------------------------------------------------------
% Represent as ABCDEF Form
%----------------------------------------------------------------

DF_model.n_s   = size(DF_model.Phi,2);
DF_model.n_eps = size(DF_model.Sigma_eta,2);
DF_model.n_y   = size(DF_model.Lambda,1);
DF_model.n_w   = size(DF_model.delta,1);
DF_model.n_e   = DF_model.n_w * DF_model.n_lags_uar;

DF_model.ABCD  = ABCD_fun_DFM(DF_model);


%% PREPARATION

%----------------------------------------------------------------
% Select Specifications
%----------------------------------------------------------------

settings.specifications = pick_var_fn(DF_model, settings);

%----------------------------------------------------------------
% Results Placeholder
%----------------------------------------------------------------

settings.est.n_methods = length(settings.est.methods_name);

partition_width = ceil(settings.simul.n_MC/num_partition); % number of MC in one partition
start_MC = (idx_partition - 1) * partition_width + 1; % start index of MC
end_MC = min(idx_partition * partition_width, settings.simul.n_MC); % end index of MC

results_irf = NaN(settings.est.n_methods,settings.est.IRF_hor,end_MC - start_MC + 1,settings.specifications.n_spec); % IRF_hor*n_MC*n_spec
results_n_lags = NaN(settings.est.n_methods,end_MC - start_MC + 1,settings.specifications.n_spec); %n_MC*n_spec

results_largest_root_svar = NaN(end_MC - start_MC + 1,settings.specifications.n_spec); % n_MC*n_spec
results_LM_stat_svar = NaN(end_MC - start_MC + 1,settings.specifications.n_spec); % n_MC*n_spec
results_LM_pvalue_svar = NaN(end_MC - start_MC + 1,settings.specifications.n_spec); % n_MC*n_spec
results_Granger_stat_svar = NaN(end_MC - start_MC + 1,settings.specifications.n_spec); % n_MC*n_spec
results_Granger_pvalue_svar = NaN(end_MC - start_MC + 1,settings.specifications.n_spec); % n_MC*n_spec
results_lambda_lp_penalize = NaN(end_MC - start_MC + 1,settings.specifications.n_spec); % n_MC*n_spec
results_weight_var_avg = NaN(2*settings.est.n_lags_max,length(settings.est.average_store_weight),...
    end_MC - start_MC + 1,settings.specifications.n_spec); % n_models*n_horizon*n_MC*n_spec
results_F_stat_svar_iv = NaN(end_MC - start_MC + 1,settings.specifications.n_spec); %n_MC*n_spec
results_F_pvalue_svar_iv = NaN(end_MC - start_MC + 1,settings.specifications.n_spec); %n_MC*n_spec


%% PRELIMINARY COMPUTATIONS: ESTIMANDS

%----------------------------------------------------------------
% Compute True IRFs in Complete Model
%----------------------------------------------------------------

[DF_model.irf, settings.est.shock_weight] = compute_irfs(DF_model,settings);

%----------------------------------------------------------------
% Compute Degree of Invertibility in Specifications
%----------------------------------------------------------------

DF_model.R0_sq = compute_invert_DFM(DF_model,settings);

%----------------------------------------------------------------
% Compute Persistency of Observables in Specifications
%----------------------------------------------------------------

[DF_model.LRV_Cov_tr_ratio, DF_model.VAR_largest_root, DF_model.frac_coef_for_large_lags] =...
    compute_persist_DFM(DF_model,settings);

%----------------------------------------------------------------
% Compute Target IRF
%----------------------------------------------------------------

switch estimand_type
    
    case 'ObsShock'
        DF_model.target_irf = DF_model.irf(settings.est.IRF_select, ...
            settings.specifications.var_select(:,settings.est.IRF_response_var_pos));
        
    case 'Recursive'
        DF_model.VAR_irf = compute_VARirfs_DFM(DF_model,settings);
        DF_model.target_irf = DF_model.VAR_irf(settings.est.IRF_select, :);
        
    case 'IV'
        DF_model.IV_irf = compute_IVirfs(DF_model,settings);
        DF_model.target_irf = DF_model.IV_irf(settings.est.IRF_select, :);

end

%----------------------------------------------------------------
% Compute Population IV Strengths
%----------------------------------------------------------------

if strcmp(estimand_type, 'IV')
    DF_model.IV_strength = compute_IVstrength_DFM(DF_model, settings);
end


%% MONTE CARLO ANALYSIS

disp('Monte Carlo simulation starts.');
disp(['dgp type: ', dgp_type]);
disp(['estimand type: ', estimand_type]);
disp(['lag type: ', num2str(lag_type)]);

parfor i_MC = 1:(end_MC - start_MC + 1)
% for i_MC = 1:settings.simul.n_MC

    if mod(i_MC, 100) == 0
        disp("Monte Carlo repetitions:")
        disp(i_MC)
    end

    %----------------------------------------------------------------
    % Generate Data
    %----------------------------------------------------------------
    
    rng(settings.simul.seed(i_MC + start_MC - 1));

    data_sim_all = generate_data(DF_model,settings);

    %----------------------------------------------------------------
    % List All Temporary Storage for i_MC in parfor
    %----------------------------------------------------------------
    
    temp_irf = NaN(settings.est.n_methods,settings.est.IRF_hor,settings.specifications.n_spec);
    temp_n_lags = NaN(settings.est.n_methods,settings.specifications.n_spec);
    
    temp_largest_root_svar = NaN(1,settings.specifications.n_spec);
    temp_LM_stat_svar = NaN(1,settings.specifications.n_spec);
    temp_LM_pvalue_svar = NaN(1,settings.specifications.n_spec);
    temp_Granger_stat_svar = NaN(1,settings.specifications.n_spec);
    temp_Granger_pvalue_svar = NaN(1,settings.specifications.n_spec);
    temp_lambda_lp_penalize = NaN(1,settings.specifications.n_spec);
    temp_weight_var_avg = NaN(2*settings.est.n_lags_max,...
        length(settings.est.average_store_weight),settings.specifications.n_spec);
    temp_F_stat_svar_iv = NaN(1,settings.specifications.n_spec);
    temp_F_pvalue_svar_iv = NaN(1,settings.specifications.n_spec);
    
    %----------------------------------------------------------------
    % Selecting Data
    %----------------------------------------------------------------

    for i_spec = 1:settings.specifications.n_spec
        
        data_sim_select = select_data_fn(data_sim_all,settings,i_spec);
    
        %----------------------------------------------------------------
        % IRF Estimation
        %----------------------------------------------------------------
        
        for i_method = 1:settings.est.n_methods
            
            switch settings.est.methods_name{i_method}

                case 'svar' % VAR
                    [temp_irf(i_method,:,i_spec),temp_n_lags(i_method,i_spec),...
                        temp_largest_root_svar(i_spec),temp_LM_stat_svar(i_spec),temp_LM_pvalue_svar(i_spec),...
                        temp_Granger_stat_svar(i_spec),temp_Granger_pvalue_svar(i_spec)]...
                        = SVAR_est(data_sim_select,settings,0);

                case 'svar_corrbias' % bias-corrected VAR
                    [temp_irf(i_method,:,i_spec),temp_n_lags(i_method,i_spec)]...
                        = SVAR_est(data_sim_select,settings,1);

                case 'bvar' % Bayesian VAR
                    [temp_irf(i_method,:,i_spec),temp_n_lags(i_method,i_spec)]...
                        = BVAR_est(data_sim_select,settings);

                case 'lp' % LP
                    [temp_irf(i_method,:,i_spec),temp_n_lags(i_method,i_spec)]...
                        = LP_est(data_sim_select,settings);

                case 'lp_penalize' % shrinkage LP
                    [temp_irf(i_method,:,i_spec),temp_n_lags(i_method,i_spec), temp_lambda_lp_penalize(i_spec)]...
                        = LP_shrink_est(data_sim_select,settings);

                case 'var_avg' % VAR model averaging
                    [temp_irf(i_method,:,i_spec),temp_n_lags(i_method,i_spec), temp_weight_var_avg(:,:,i_spec)]...
                        = VAR_avg_est(data_sim_select,settings);

                case 'svar_iv' % SVAR-IV       
                    [temp_irf(i_method,:,i_spec),temp_n_lags(i_method,i_spec),...
                        temp_F_stat_svar_iv(i_spec),temp_F_pvalue_svar_iv(i_spec)]...
                        = SVAR_IV_est(data_sim_select,settings);
                
            end
            
        end
        
    end
    
    %----------------------------------------------------------------
    % Move Results to Permanent Storage in parfor
    %----------------------------------------------------------------
    
    results_irf(:,:,i_MC,:) = temp_irf;
    results_n_lags(:,i_MC,:) = temp_n_lags;
    
    results_largest_root_svar(i_MC,:) = temp_largest_root_svar;
    results_LM_stat_svar(i_MC,:) = temp_LM_stat_svar;
    results_LM_pvalue_svar(i_MC,:) = temp_LM_pvalue_svar;
    results_Granger_stat_svar(i_MC,:) = temp_Granger_stat_svar;
    results_Granger_pvalue_svar(i_MC,:) = temp_Granger_pvalue_svar;
    results_lambda_lp_penalize(i_MC,:) = temp_lambda_lp_penalize;
    results_weight_var_avg(:,:,i_MC,:) = temp_weight_var_avg;
    results_F_stat_svar_iv(i_MC,:) = temp_F_stat_svar_iv;
    results_F_pvalue_svar_iv(i_MC,:) = temp_F_pvalue_svar_iv;

end

% clear temporary storage
clear temp_* i_MC i_spec data_sim_all data_sim_select i_method ...
    partition_width start_MC end_MC


%% SUMMARIZE RESULTS

%----------------------------------------------------------------
% Wrap up Results, Pick out Target IRF
%----------------------------------------------------------------

% wrap up results from parallel loop

for i_method = 1:settings.est.n_methods
    
    thisMethod = settings.est.methods_name{i_method};
    results.irf.(thisMethod) = permute(results_irf(i_method,settings.est.IRF_select,:,:), [2 3 4 1]);
    results.n_lags.(thisMethod) = permute(results_n_lags(i_method,:,:), [2 3 1]);
    
end

if any(strcmp(settings.est.methods_name, 'svar'))    
    results.largest_root.svar = results_largest_root_svar;
    results.LM_stat.svar = results_LM_stat_svar;
    results.LM_pvalue.svar = results_LM_pvalue_svar;
    if strcmp(estimand_type, 'IV')
        results.Granger_stat.svar = results_Granger_stat_svar;
        results.Granger_pvalue.svar = results_Granger_pvalue_svar;
    end
end

if any(strcmp(settings.est.methods_name, 'lp_penalize'))    
    results.lambda.lp_penalize = results_lambda_lp_penalize;
end

if any(strcmp(settings.est.methods_name, 'var_avg'))
    results.weight.var_avg = results_weight_var_avg;
end

if any(strcmp(settings.est.methods_name, 'svar_iv'))
    results.F_stat.svar_iv = results_F_stat_svar_iv;
    results.F_pvalue.svar_iv = results_F_pvalue_svar_iv;
end

clear results_* i_method thisMethod

%----------------------------------------------------------------
% Compute Mean-Squared Errors, Bias-Squared, Variance
%----------------------------------------------------------------

% compute MSE, Bias2, VCE for each horizon and each specification

[results.MSE, results.BIAS2, results.VCE] = irf_perform_summary(results.irf, DF_model.target_irf, settings);

%----------------------------------------------------------------
% Export Results
%----------------------------------------------------------------

% save results or one partition of results

mkdir(save_folder); 
if num_partition == 1
    save(fullfile(save_folder, strcat('DFM_', dgp_type, '_', estimand_type)), ...
    'DFM_estimate','DF_model','settings','results','-v7.3'); % save results
else
    save(fullfile(save_folder, strcat('DFM_', dgp_type, '_', estimand_type, '_', num2str(idx_partition))), ...
    'DFM_estimate','DF_model','settings','results','-v7.3'); % save one partition of results
    partition_log_file = fopen(fullfile(save_folder, 'partition_log.txt'),'a+'); % report idx_partition in log file once finish exporting
    fprintf(partition_log_file, '%d\n', idx_partition);
    fclose(partition_log_file);
end
clear idx_partition

% combine results for multiple partitions

if num_partition > 1
    
    % detect if all the partitions have finished
    
    partition_log_file = fopen(fullfile(save_folder, 'partition_log.txt'),'r');
    partition_finished_list = fscanf(partition_log_file, '%d'); % read in all the indices for finished partitions
    fclose(partition_log_file);
    
    % start to combine when all the partitions have finished
    
    if length(partition_finished_list) == num_partition       
        clear results
        results = combine_results(save_folder, strcat('DFM_', dgp_type, '_', estimand_type), num_partition);
        
        % recompute MSE, Bias2, VCE
        [results.MSE, results.BIAS2, results.VCE] = irf_perform_summary(results.irf, DF_model.target_irf, settings);
        
        % save combined results
        save(fullfile(save_folder, strcat('DFM_', dgp_type, '_', estimand_type)), ...
            'DFM_estimate','DF_model','settings','results','-v7.3');
    end
    
end

delete(poolobj);
clear save_folder save_suff num_partition partition_* poolobj
toc;


%% PLOT RESULTS

%----------------------------------------------------------------
% Plot IRFs for Checking
%----------------------------------------------------------------

% for i_method = 1:settings.est.n_methods
%     
%     thisMethod = settings.est.methods_name{i_method};
%     figure(i_method)
%     plot(settings.est.IRF_select, DF_model.target_irf(:,settings.specifications.plot_indx),'Linewidth',5)
%     hold on
%     for i = 1:settings.simul.n_MC
%         eval(['plot(settings.est.IRF_select, results.irf.' thisMethod '(:,i,settings.specifications.plot_indx))']);
%         hold on
%     end
%     title(replace(thisMethod,'_',' '))
%     hold off
% 
% end
% 
% clear i_method thisMethod i