%% Conditional VAR (Bayesian Estimation via MCMC)
% This script estimates a conditional VAR model using Bayesian MCMC methods,
% with emphasis on tail and conditional risk transmission.

%% 1. Housekeeping
clear; close all; clc; delete *asv;

% VAR lag order
p = 2;

% MCMC settings
Nsim  = 1e5;          % Total number of simulations
Nburn = 0.2 * Nsim;  % Burn-in period

%% 2. Load and Transform Data
% Import raw macro-financial data
dataraw = importdata('PATH_TO_DATA/data_new.xlsx');

% Transform GDP series to growth rate (log-difference, percentage)
datapart = 100 * log(dataraw(2:end,1) ./ dataraw(1:end-1,1));
datapart = [-6.275; datapart];   % Initial value adjustment

% Replace original GDP series
dataraw(:,1) = datapart;
data = dataraw;

% Initial conditions (first four observations)
Y0 = data(1:4, :);

% Effective sample
Y = data(5:end, :);
[T, n] = size(Y);

% Stack observations for SUR representation
y = reshape(Y', T * n, 1);

%% 3. Construct Regressor Matrix
% Number of coefficients per equation (including intercept)
k = 7;

% Create lagged regressor matrix
tmpY = [Y0(end-1:end, :); Y];
X_tilde = zeros(T, 6);   % Excluding intercept initially

for i = 1:2
    X_tilde(:, (i-1)*n+1 : i*n) = tmpY(2-i+1:end-i, :);
end

% Add intercept
X_tilde = [ones(T,1), X_tilde];

% Convert to SUR form
x = SURform2(X_tilde, n);

%% 4. Prior Specification
% Prior for covariance matrix (Inverse-Wishart)
nu0 = n + 3;
S0  = eye(n);

% Prior mean for VAR coefficients
beta0 = zeros(n * k, 1);

% Prior precision matrix
% Intercept variance set lower than slope coefficients
tmp = ones(k * n, 1);
tmp(1:2*3+1:k*n) = 1/10;
ivbeta = sparse(1:k*n, 1:k*n, tmp);

%% 5. Storage for Posterior Draws
store_sig  = zeros(Nsim, n, n);
store_beta = zeros(Nsim, n * k);

%% 6. Initialization (MLE Starting Values)
% Initial OLS estimates
betta = (x' * x) \ (x' * y);

% Residuals and covariance initialization
e    = reshape(y - x * betta, n, T);
Sig  = e * e' / T;
iSig = Sig \ speye(n);

%% 7. MCMC Sampling
for ii = 1 : (Nsim + Nburn)

    % Step 1: Sample beta | Sigma, Y  ~ Normal
    Xisig   = x' * kron(speye(T), iSig);
    XisigX  = Xisig * x;
    Dbeta   = ivbeta + XisigX;
    XisigY  = Xisig * y;

    betahat = Dbeta \ (ivbeta * beta0 + XisigY);
    betta   = betahat + chol(Dbeta, 'lower') * randn(21,1);

    % Step 2: Sample Sigma | beta, Y  ~ Inverse-Wishart
    e    = reshape(y - x * betta, n, T);
    Sig  = iwishrnd(S0 + e * e', nu0 + T);
    iSig = Sig \ speye(n);

    % Progress indicator
    if mod(ii, 1000) == 0
        disp([num2str(ii), ' loops completed...']);
    end

    % Store posterior draws after burn-in
    if ii > Nburn
        jj = ii - Nburn;
        store_beta(jj, :)   = betta';
        store_sig(jj, :, :) = Sig;
    end
end

%% 8. Posterior Inference
% Posterior moments
beta_mean = mean(store_beta);
beta_var  = var(store_beta);

Sig_mean = squeeze(mean(store_sig));
Sig_var  = squeeze(var(store_sig));

fprintf('Posterior mean of beta:\n');
disp(beta_mean);

fprintf('Posterior variance of beta:\n');
disp(beta_var);

fprintf('-----------------------------------------------\n');

fprintf('Posterior mean of Sigma:\n');
disp(Sig_mean);

fprintf('Posterior variance of Sigma:\n');
disp(Sig_var);

fprintf('-----------------------------------------------\n');

% Credible intervals (95%)
beta_CI    = quantile(store_beta, [0.025 0.975]);
betavar_CI = quantile(beta_var,  [0.025 0.975]);
Sig_CI     = quantile(store_sig, [0.025 0.975]);
sigvar_CI  = quantile(Sig_var,   [0.025 0.975]);
