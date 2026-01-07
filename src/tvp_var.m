%% Time-Varying Parameter VAR (TVP-VAR)
% Bayesian estimation of a TVP-VAR model using MCMC,
% with time-specific impulse response analysis.

%% 1. Housekeeping
clear; close all; clc; delete *asv;

% MCMC settings
Nsim  = 1e5;           % Total number of simulations
Nburn = 0.2 * Nsim;   % Burn-in period

% Model dimensions
nk   = 21;    % Number of coefficients per time period
n_hz = 10;    % Impulse response horizon

% Selected time points for IRF comparison
tt1 = 33;     % Early sample (e.g. 1975Q1)
tt2 = 65;     % Crisis period (e.g. 2008Q1)

%% 2. Load and Transform Data
% Import raw macro-financial data
dataraw = importdata('PATH_TO_DATA/data_new.xlsx');

% Transform GDP into growth rate (log-difference, percentage)
datapart = 100 * log(dataraw(2:end,1) ./ dataraw(1:end-1,1));
datapart = [-6.275; datapart];

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
% Lagged regressors (two lags)
X_tilde = zeros(T, n * 2);
tmpY = [Y0(end-2+1:end, :); Y];

for i = 1:2
    X_tilde(:, (i-1)*n+1 : i*n) = tmpY(2-i+1:end-i, :);
end

% SUR-form design matrix with intercept
X = SURform([ones(n*T,1), kron(X_tilde, ones(n,1))]);

% State evolution matrix for time-varying coefficients
Hb = speye(T*nk, T*nk) - ...
     sparse(nk+1:T*nk, 1:(T-1)*nk, ones(1,(T-1)*nk), T*nk, T*nk);

%% 4. Prior Specification
% Prior for covariance matrix Sigma
nuSig0 = n + 3;
SSig0  = eye(n);

% Prior mean for initial VAR coefficients
mubeta0 = zeros(nk,1);

% Prior precision for beta (intercept treated separately)
tmp = ones(nk,1);
tmp(1:3*2+1:nk) = 1/10;
iVbeta = sparse(1:nk, 1:nk, tmp);

% Prior for state innovation variances Q
nuQ0 = 3 * ones(nk,1);
QS0  = 0.01^2 * ones(nk,1);
QS0(1:3*2+1:nk) = 0.1^2;

%% 5. Storage Initialization
store_beta  = zeros(Nsim, T*nk);
store_Q     = zeros(Nsim, nk);
store_Sigma = zeros(Nsim, n, n);

% Storage for impulse responses
store_yIR11 = zeros(n_hz, n);
store_yIR12 = zeros(n_hz, n);
store_yIR13 = zeros(n_hz, n);
store_yIR21 = zeros(n_hz, n);
store_yIR22 = zeros(n_hz, n);
store_yIR23 = zeros(n_hz, n);

%% 6. Initialization (OLS Starting Values)
Z = SURform2([ones(T,1), X_tilde], n);
beta0 = (Z' * Z) \ (Z' * y);

e    = reshape(y - Z * beta0, n, T);
Sig  = e * e' / T;
iSig = Sig \ speye(n);

Q = 0.1 * ones(nk,1);

%% 7. MCMC Sampling
for ii = 1 : (Nsim + Nburn)

    % 1. Sample time-varying beta
    Xisig  = X' * kron(speye(T), iSig);
    XisigX = Xisig * X;

    iQ = sparse(1:T*nk, 1:T*nk, repmat(1 ./ Q, T, 1));
    Dbeta = iQ + XisigX;

    XisigY = Xisig * y;
    tildealpha = kron(ones(T,1), beta0);

    betahat = Dbeta \ (iQ * tildealpha + XisigY);
    betta   = betahat + chol(Dbeta, 'lower') * randn(nk*T,1);

    % 2. Sample Sigma
    e    = reshape(y - X * betta, n, T);
    Sig  = iwishrnd(SSig0 + e * e', nuSig0 + T);
    iSig = Sig \ speye(n);

    % 3. Sample Q (state innovation variances)
    nuu = nuQ0 + T/2;
    u   = reshape(betta - [beta0; betta(1:end-nk)], nk, T);
    Q   = 1 ./ gamrnd(nuu, 1 ./ (sum(u.^2,2)/2 + QS0));

    % 4. Sample initial beta0
    ivQ = sparse(1:nk, 1:nk, 1 ./ Q);
    Dbeta0 = ivQ + iVbeta;

    beta0hat = Dbeta0 \ (iVbeta * mubeta0 + ivQ * betta(1:21));
    beta0    = beta0hat + chol(Dbeta0, 'lower')' * randn(21,1);

    % Store posterior draws and impulse responses
    if ii > Nburn
        jj = ii - Nburn;

        store_beta(jj,:)    = betta';
        store_Q(jj,:)       = Q';
        store_Sigma(jj,:,:) = Sig;

        % Structural shocks (Cholesky identification)
        CSig = chol(Sig, 'lower');
        shock1 = [1;0;0] / CSig(1,1);
        shock2 = [0;1;0] / CSig(2,2);
        shock3 = [0;0;1] / CSig(3,3);

        tempbeta = reshape(betta, nk, T);

        % Impulse responses at selected time points
        imp1_obj1 = construct_IR(tempbeta(:,tt1), Sig, n_hz, shock1);
        imp2_obj1 = construct_IR(tempbeta(:,tt1), Sig, n_hz, shock2);
        imp3_obj1 = construct_IR(tempbeta(:,tt1), Sig, n_hz, shock3);

        imp1_obj2 = construct_IR(tempbeta(:,tt2), Sig, n_hz, shock1);
        imp2_obj2 = construct_IR(tempbeta(:,tt2), Sig, n_hz, shock2);
        imp3_obj2 = construct_IR(tempbeta(:,tt2), Sig, n_hz, shock3);

        store_yIR11 = store_yIR11 + imp1_obj1;
        store_yIR12 = store_yIR12 + imp2_obj1;
        store_yIR13 = store_yIR13 + imp3_obj1;

        store_yIR21 = store_yIR21 + imp1_obj2;
        store_yIR22 = store_yIR22 + imp2_obj2;
        store_yIR23 = store_yIR23 + imp3_obj2;
    end
end

%% 8. Posterior Mean Impulse Responses
yIR11_mean = store_yIR11 / Nsim;
yIR12_mean = store_yIR12 / Nsim;
yIR13_mean = store_yIR13 / Nsim;

yIR21_mean = store_yIR21 / Nsim;
yIR22_mean = store_yIR22 / Nsim;
yIR23_mean = store_yIR23 / Nsim;

%% 9. Visualization: Aggregate IRFs
figure('color','w')
subplot(2,3,1); plot(yIR11_mean); box off; xlim([1,n_hz]);
title('GDP shock impulse (1975Q1)');

subplot(2,3,2); plot(yIR12_mean); box off; xlim([1,n_hz]);
title('FFR shock impulse (1975Q1)');

subplot(2,3,3); plot(yIR13_mean); box off; xlim([1,n_hz]);
title('Inflation shock impulse (1975Q1)');

subplot(2,3,4); plot(yIR21_mean); box off; xlim([1,n_hz]);
title('GDP shock impulse (2008Q1)');

subplot(2,3,5); plot(yIR22_mean); box off; xlim([1,n_hz]);
title('FFR shock impulse (2008Q1)');

subplot(2,3,6); plot(yIR23_mean); box off; xlim([1,n_hz]);
title('Inflation shock impulse (2008Q1)');

%% 10. Posterior Summary Statistics
beta_mean = mean(store_beta);
Sig_mean  = squeeze(mean(store_Sigma));
Q_mean    = mean(store_Q);

beta_var = var(store_beta);
Sig_var  = squeeze(var(store_Sigma));
Q_var    = var(store_Q);

% Credible intervals
beta_CI = quantile(store_beta, [0.025 0.975]);
Sig_CI  = quantile(store_Sigma, [0.025 0.975]);
Q_CI    = quantile(store_Q,     [0.025 0.975]);
