%%% clear everything
close all
clear
clc

%%% check if CVX is installed
try
    run('cvx_setup');
catch err
    error('CVX problem.');
end

%%% setup
% number of sensors
m = 100;
% size of the unknown
n = 20;

% measurement matrix of the sensor network, tight
A = randn(n, m);
[U, ~, V] = svd(A, 'econ');
A = sqrt(m)*U*V';

% best achievable performance in terms of MSE
minMSE = trace_inv(A*A');
% target performance
rho = 2;

%%% select a small number of sensors that achieve the prescribed MSE
[z, time1] = select_mse(A, rho*minMSE);

%%% select sensors uniformly over T time instances
% number of time instances
T = 10;
% regularization parameter
lambda = 1000;
% call
[Z1, time2] = select_mse_implicit_energy(A, T, rho*minMSE, lambda);

%%% select sensors with explicit power constraints
% number of time instances
T = 10;
% regularization parameter
lambda = 1000;
% sensing costs
s = norms(A).^2;
% network costs
load('theC.mat');
% reference energy levels, if available
e0 = zeros(m,1);
% call
[Z2, e, time3] = select_mse_explicit_energy(A, T, rho*minMSE, lambda, s, C, e0);
