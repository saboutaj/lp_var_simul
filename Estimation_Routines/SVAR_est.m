function [IRF,nlags,largest_root,LM_stat] = SVAR_est(data_sim,settings);

% preparations

IRF_hor    = settings.est.IRF_hor;
n_lags_max = settings.est.n_lags_max;
est_n_lag  = settings.est.est_n_lag;
est_n_lag_BIC  = settings.est.est_n_lag_BIC;
n_lags_fix = settings.est.n_lags_fix;

res_autocorr_nlags = 4; % check autocorr of residuals up to this order

response_pos = settings.est.IRF_response_var_pos;

with_shock = settings.est.with_shock;
recursive_shock = settings.est.recursive_shock;
with_IV = settings.est.with_IV;

if recursive_shock == 1
    recursive_shock_pos = settings.est.recursive_shock_pos;
end
if with_IV == 1
    IV_est_normalize_var_pos = settings.est.IV_est_normalize_var_pos;
end

% collect data

if with_shock == 1
    Y = [data_sim.data_shock,data_sim.data_y];
    responseV = response_pos + 1;
    recursiveShock = 1;
    normalizeV = recursiveShock;
elseif with_IV == 1
    Y = [data_sim.data_z,data_sim.data_y];
    responseV = response_pos + 1;
    recursiveShock = 1;
    normalizeV = IV_est_normalize_var_pos + 1;
else
    Y = data_sim.data_y;
    responseV = response_pos;
    recursiveShock = recursive_shock_pos;
    normalizeV = recursiveShock;
end

% set lag length

if est_n_lag == 0
    nlags = n_lags_fix;
elseif est_n_lag_BIC == 1
    [BIC,~] = IC_VAR(Y,n_lags_max);
    [~,nlags] = min(BIC);
else
    [~,AIC] = IC_VAR(Y,n_lags_max);
    [~,nlags] = min(AIC);
end

% estimate VAR

[~,By,Sigma,~,Res] = VAR(Y,nlags);
G = chol(Sigma, 'lower');
ShockVector = G(:,recursiveShock);

% estimate IRF

IRF = IRF_SVAR(By,ShockVector,IRF_hor - 1);
IRF = IRF(responseV,:) / IRF(normalizeV,1);
IRF = IRF';

% estimate largest root in VAR

nv = size(Y,2);
companion_form = diag(ones(1, nv*(nlags-1)), -nv);
companion_form(1:nv,:) = reshape(By,[nv, nv*nlags]);
largest_root = max(abs(eig(companion_form)));

% estimate LM-stat to examine VAR(p) fit

nT = size(Y,1);
Y_lag = lagmatrix(Y,1:nlags); % lagged Y as explanatory variables
Y_lag = Y_lag((nlags+res_autocorr_nlags+1):end,:);
Res_lag = lagmatrix(Res,1:res_autocorr_nlags); % lagged residual
Res_lag = Res_lag((res_autocorr_nlags+1):end, :);
X_auxiliary = [ones(nT-nlags-res_autocorr_nlags,1), Y_lag, Res_lag];
Res_current = Res((res_autocorr_nlags+1):end, :); % current residual
[~,~,~,Res_aux] = LS(Res_current, X_auxiliary);
Sigma_original = cov(Res_current,1);
Sigma_auxiliary = cov(Res_aux,1);
LM_stat = (nT-nlags-res_autocorr_nlags - nv*nlags - 1 - nv*res_autocorr_nlags - 0.5) *...
    log(det(Sigma_original) / det(Sigma_auxiliary));

end