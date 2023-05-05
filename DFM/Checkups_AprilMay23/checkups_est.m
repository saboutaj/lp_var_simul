%%rerun the estimation on a long-run simulated dataset with parameters estimated from the
%%SW dataset - checks if we recover the parameters
clear all;

addpath '/home/sa6128/Desktop/lp_var_simul/DFM/Subroutines/SW_DFM_Estimation/matlab'
addpath '/home/sa6128/Desktop/lp_var_simul/DFM/Subroutines/SW_DFM_Estimation/m_utilities'
addpath '/home/sa6128/Desktop/lp_var_simul/Estimation_Routines'
addpath '/home/sa6128/Desktop/lp_var_simul/DFM/Settings'

n_y = 207; % number of macro-variables
reorder    = [1:76, 87:94, 77:86, 95:171, 181:195, 172:180, 196:207]; % index to reorder data to match variable list in Stock-Watson (2016)
levels     = 1; % =1: variables in levels, 0=: differenced
n_fac      = 6; % number of factors
n_lags_fac = 2; % lag order of factors
n_lags_uar = 2; % lag order of measurement error
spec_id = 1; % seed for random draws of specifications (= DGPs from encompassing model)
lag_type = 2; 
mode_type = 1; % robustness check mode:
               % 1 (baseline), 2 (cumulative IRF), 
               % 3 (persistent DGP), 4 (persistent DGP with MN prior), 
               % 5 (small sample), 6 (salient series)
estim_diagn = 1; % =1: show DFM estimation diagnostics

DFM = DFM_est(n_fac ,n_lags_fac ,n_lags_uar , reorder, levels); %saving the estimated parameters

%% simulate the estimated DGP
% Settings
rng('default');
shared;
%G;
%ObsShock;
%check_mode;
DF_model.levels = levels;
settings.DFM = DF_model;
settings.simul.T = 10000;
settings.simul.T_burn = settings.simul.T/100;
repeat_sim = 30;

%Initiate Simulation Data
sim_SW = zeros(settings.simul.T,n_y,repeat_sim);%time series simulation data
DFM_simdata(repeat_sim+1) = struct();%estimation data
fields = fieldnames(DFM);
fields = fields(1:6);%keep only parameter fields

%initiate the fields for parameters and the bias and Std matrices
for j=1:length(fields)
    field = fields{j};
    DFM_simdata(1).(field)= DFM.(field); %first line of estimation data are the real coefficients
    bias.(field) = -(repeat_sim)*DFM_simdata(1).(field);
    Std.(field) = 0;
end

%Need to parallelize
datain(repeat_sim) = struct();%estimation data


%%
tic
%run the simulation
parfor i = 1:repeat_sim
    [sim_SW(:,:,i),shocks,DFM_simdata(i+1).realFac] = generate_data_SW(DFM,settings);
% rerun the estimation on the simulation
    datain(i).bpdata = sim_SW(:,:,i);      
    datain(i).bpinclcode = DFM.bpinclcode;
    datain(i).bpnamevec = DFM.bpnamevec;
    datain(i).bplabvec_long = DFM.bplabvec_long;
    datain(i).bplabvec_short = DFM.bplabvec_short;
    datain(i).bptcodevec = DFM.bptcodevec;
    %adapt calvec
    [dnobs_sim,calvec_sim,calds_sim] = calendar_make([1959 1],[1959+(settings.simul.T/4)-1 4],4);
    datain(i).calvec = calvec_sim;
    save_DFM = DFM_est_sim(n_fac ,n_lags_fac ,n_lags_uar , reorder, levels,datain(i),settings.simul.T);
    for j=1:length(fields)
        field = fields{j};
        DFM_simdata(i+1).(field)= save_DFM.(field);
    end %saving the estimated parameters
    DFM_simdata(i+1).estimFac= save_DFM.fac;
end
toc
%% Adapt parameters by computing the rotation
Y = [DFM_simdata(2).estimFac];
X = [DFM_simdata(2).realFac];
for i=3:repeat_sim+1
    Y = vertcat(Y,[DFM_simdata(i).estimFac]);
    X = vertcat(X,[DFM_simdata(i).realFac]);
end
if levels==1
    timetrend = repmat([1:settings.simul.T],[1,repeat_sim])';
    X = [timetrend,X];
end
[beta,sigma] = mvregress(X,Y);
if levels==1
    H = beta(2:7,:)';
else 
    H = beta(:,:)';
end
invH = inv(H);
for i=1:repeat_sim
    DFM_simdata(i+1).Lambda = DFM_simdata(i+1).Lambda*H;
    DFM_simdata(i+1).Phi(1:6,:) = horzcat(invH*DFM_simdata(i+1).Phi(1:6,1:6)*H,invH*DFM_simdata(i+1).Phi(1:6,7:12)*H);
    DFM_simdata(i+1).Sigma_eta = invH*DFM_simdata(i+1).Sigma_eta*(invH)';
end
%% Compute bias and standard deviations
for j=1:length(fields)
    field = fields{j};
    for i=1:repeat_sim
        bias.(field) = bias.(field)+DFM_simdata(i+1).(field) ;
        Std.(field) = Std.(field)+(DFM_simdata(i+1).(field)).^2 ;
    end
    bias.(field) = bias.(field)/repeat_sim;
    Std.(field) = sqrt((Std.(field)/repeat_sim)-(bias.(field)+DFM_simdata(1).(field)).^2);
    writematrix([bias.(field)],"Bias_level="+string(levels)+".xls",'Sheet',j)
    writematrix([Std.(field)],"Std_level="+string(levels)+".xls",'Sheet',j)
    writematrix([DFM_simdata(1).(field)],"Parameter_simul_level="+string(levels)+".xls",'Sheet',j)
end
bias.LambdaPct = 100*bias.Lambda./DFM.Lambda;
bias.sigmavPct = 100*bias.sigma_v./DFM.sigma_v;
bias.phiPct = 100*bias.Phi./DFM.Phi;
bias.deltaPct = 100*bias.delta./DFM.delta;
bias.sigmaetaPct = 100*bias.Sigma_eta./DFM.Sigma_eta;
bias.LambdaAv = mean([bias.Lambda]);
bias.sigmaAv = mean([bias.sigma_v]);
bias.deltaAv = mean([bias.delta]);
bias.LambdaMax = max(abs([bias.Lambda]));
bias.sigmaMax = max(abs([bias.sigma_v]));
bias.deltaMax = max(abs([bias.delta]));
writematrix([bias.LambdaPct],"Bias_level="+string(levels)+".xls",'Sheet',j+1)
writematrix([bias.phiPct],"Bias_level="+string(levels)+".xls",'Sheet',j+2)
writematrix([bias.sigmaetaPct],"Bias_level="+string(levels)+".xls",'Sheet',j+3)
writematrix([bias.sigmavPct],"Bias_level="+string(levels)+".xls",'Sheet',j+4)
writematrix([bias.deltaPct],"Bias_level="+string(levels)+".xls",'Sheet',j+5)
%writematrix([bias.deltaAv ; bias.deltaMax],"Bias2.xls",'Sheet',j+4)