%% Load Data
% Import macro-financial time series data from Excel file
data = importdata('PATH_TO_DATA/data_new.xlsx');

% Construct quarterly time index
tt = (1967:0.25:2022.5)';

% Plot raw data
figure('color','w');
plot(tt, [data], 'linewidth', 3);
xlim([1967 2023]);

%% Stationarity Testing: Augmented Dickey-Fuller (ADF)
% Perform ADF tests on level data (monthly series)
ans1 = adftest(data(:,1));  % GDP
ans2 = adftest(data(:,2));  % Federal Funds Rate (FFR)
ans3 = adftest(data(:,3));  % Inflation

% First differences for non-stationary series
dgdp = diff(data(:,1));
dr   = diff(data(:,2));

% ADF tests on differenced series
ansi = adftest(dgdp);
ansr = adftest(dr);

%% Autocorrelation and Partial Autocorrelation Analysis
% Compute autocorrelation functions
gdp_ar = autocorr(data(:,1));
r_ar   = autocorr(data(:,2));
inf_ar = autocorr(data(:,3));

% Compute partial autocorrelation functions
gdp_pr = parcorr(data(:,1));
r_pr   = parcorr(data(:,2));
inf_pr = parcorr(data(:,3));

% Plot ACF and PACF for all variables
subplot(2,3,1)
stem(gdp_ar);
title('GDP ACF');

subplot(2,3,2)
stem(r_ar);
title('FFR ACF');

subplot(2,3,3)
stem(inf_ar);
title('Inflation ACF');

subplot(2,3,4)
stem(gdp_pr);
title('GDP PACF');

subplot(2,3,5)
stem(r_pr);
title('FFR PACF');

subplot(2,3,6)
stem(inf_pr);
title('Inflation PACF');

%% Data Transformation
% Reload raw data for transformation
dataraw = importdata('PATH_TO_DATA/data_new.xlsx');

% Convert GDP to growth rate (log difference, annualized percentage)
datapart = 100 * log(dataraw(2:end,1) ./ dataraw(1:end-1,1));
datapart = [-6.275; datapart];  % Initial value adjustment

% Replace original GDP series with transformed series
dataraw(:,1) = datapart;
data = dataraw;

% Plot transformed data
figure('color','w');
plot(tt, [data], 'linewidth', 3);
xlim([1967 2023]);
