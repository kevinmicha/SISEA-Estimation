% Detection and Estimation - Simulations corresponding to exercises 6 to 8
% Master SISEA
% By Kevin Michalewicz
clear all
close all
clc 

N_vec = [10 50 100 500 1000];
A = 1;
sigma_vec = [0.25 0.5 1];
freq0 = 0.2; % between 0 and 0.5
phi = 1;
K = 20; % how many estimators for each (sigma, N)
phi_ML = zeros(1,K);
A_ML = zeros(size(phi_ML));
A_squared_error = zeros(length(sigma_vec),length(N_vec));
phi_squared_error = zeros(size(A_squared_error));
A_bias = zeros(length(sigma_vec),length(N_vec));
phi_bias = zeros(size(A_bias));

N_idx = 4;
sigma_idx = 2;
A_fisher = reshape((1./(2*sigma_vec.^2))'*N_vec, [], 1);
[~,fisher_order] = sort(A_fisher); % getting sortIdx (increasing order)
phi_fisher = A_fisher*A^2;

for i=1:length(sigma_vec)
    sigma = sigma_vec(i);
    for j=1:length(N_vec)
        for k=1:K
            N = N_vec(j);
            n = 1:N;
            S_N = sin(2*pi*freq0*n);
            C_N = cos(2*pi*freq0*n);
            X_N = A*cos(2*pi*freq0*n+phi) + randn(1,N)*sigma;
            phi_ML(1,k) = atan(-X_N*S_N'/(X_N*C_N'));
            A_ML(1,k) = 2*mean(X_N.*cos(2*pi*freq0*n+phi_ML(1,k)));
        end
        A_squared_error(i,j) = mean((A - A_ML).^2); 
        A_bias(i,j) = A - mean(A_ML);
        phi_squared_error(i,j) = mean((phi - phi_ML).^2);
        phi_bias(i,j) = phi - mean(phi_ML);
    end
end

% Plotting part of exercise 7

figure(1)
plot(N_vec, A_squared_error(sigma_idx,:), '-o')
hold on
plot(N_vec, A_bias(sigma_idx,:), '-o')
legend('Squared error', 'Bias');
title("Amplitude estimation for \sigma = " + sigma_vec(sigma_idx));
xlabel("Number of points");
ylabel("Error/bias value");

figure(2)
plot(sigma_vec, A_squared_error(:,N_idx), '-o')
hold on
plot(sigma_vec, A_bias(:,N_idx), '-o')
legend('Squared error', 'Bias');
title("Amplitude estimation for N = " + N_vec(N_idx));
xlabel("Standard deviation");
ylabel("Error/bias value");

figure(3)
plot(N_vec, phi_squared_error(sigma_idx,:), '-o')
hold on
plot(N_vec, phi_bias(sigma_idx,:), '-o')
legend('Squared error', 'Bias');
title("Phase estimation for \sigma = " + sigma_vec(sigma_idx));
xlabel("Number of points");
ylabel("Error/bias value");

figure(4)
plot(sigma_vec, phi_squared_error(:,N_idx), '-o')
hold on
plot(sigma_vec, phi_bias(:,N_idx), '-o')
legend('Squared error', 'Bias');
title("Phase estimation for N = " + N_vec(N_idx));
xlabel("Standard deviation");
ylabel("Error/bias value");

% Exercise 8

% Flattening the quadratic errors and bias
A_squared_error = reshape(A_squared_error,[],1);
A_bias = reshape(A_bias,[],1);
phi_squared_error = reshape(phi_squared_error,[],1);
phi_bias = reshape(phi_bias,[],1);

figure(5)
plot(A_fisher(fisher_order), A_squared_error(fisher_order), '-o')
hold on
plot(A_fisher(fisher_order), A_bias(fisher_order), '-o')
legend('Squared error', 'Bias');
title("Amplitude estimation");
xlabel("Fisher information");
ylabel("Error/bias value");

figure(6)
plot(phi_fisher(fisher_order), phi_squared_error(fisher_order), '-o')
hold on
plot(phi_fisher(fisher_order), phi_bias(fisher_order), '-o')
legend('Squared error', 'Bias');
title("Phase estimation");
xlabel("Fisher information");
ylabel("Error/bias value");
