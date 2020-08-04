clear all; clc;

% Geometric Sequence Decomposition with k-simplxes Transform.
% Author: Woong-Hee Lee, Jong-Ho Lee, and Ki Won Sung
% Date: 2020-08-02

k = 6; % Number of geometric sequence
IT = randn(1,k) + 1i*randn(1,k); % Complex-valued initial terms
CR = randn(1,k) + 1i*randn(1,k); % Complex-valued common ratios
P = 20; % Length of geometric sequence (it must be larger than 2k)

s = IT*transpose(CR).^([0:P-1]); % Superposition of geometric sequences

ek = GSDST_kEstimator(s);             % Estimation of k
[eIT, eCR] = GSDST_solver(1:ek,s,ek); % Estimation of IT and CR
% In eIT and eCR, the values located at the same index form their own
% geometric sequence like original IT and CR.

% Reconstruct s
sr = eIT * transpose(eCR).^([0:P-1]);

% Check reconstruction error
norm(s - sr) / norm(s)