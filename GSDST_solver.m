function [eIT, eCR] = GSDST_solver(sampleVec,observation, ek)

% Geometric Sequence Decomposition with k-simplxes Transform.
% Author: Woong-Hee Lee, Jong-Ho Lee, and Ki Won Sung
% Date: 2020-08-02
% 
% eIT: Estimated initial terms
% eCR: Estimated common ratios
% sampleVec: Indices for 1st vertex whose length is ek
% observation: Superposition of ek number of geometric sequences
% ek: Estimated number of geometric sequences

% Example (ek=k=3)
% k = 3;
% IT = randn(1,k); % Initial term
% CR = randn(1,k); % Common ratio
% P = 12; % Number of geometric sequence
% s = IT*transpose(CR).^([0:P-1]); % Superposition of k geometric sequences
% ek = GSDST_kEstimator(s);             % Estimation of k
% [eIT, eCR] = GSDST_solver([1 2 3], s, ek); % Most normal case
% [eIT, eCR] = GSDST_solver([1 2 9], s, ek); % Non-consecutive sampling

    % Make search space
    vMat = [];
    for k=1:ek+1
        vMat = [vMat transpose(observation(sampleVec+k-1))];
    end
    
    % Make simplexes and polynomials
    candiV = combinator(ek+1,ek,'c');
    for k=1:ek+1
        selecIndex = candiV(k,:);
        pMat = [];
        for kk=1:length(selecIndex)
            pMat = [pMat vMat(:,selecIndex(kk))];
        end
        if mod(k,2)~=0
            p(k) = det(pMat);
        else
            p(k) = -det(pMat);
        end
    end
    
    % Find roots of p, i.e. find eCR
    sol = transpose(roots(p));
    % One can use another function instead of the roots function to increase accuracy.
        
    % Find eIT
    T = zeros(ek,ek);
    for k=1:ek
        for kk=1:ek
       T(k,kk) = sol(k).^(sampleVec(kk)-1); 
        end
    end
    eIT = observation(sampleVec)*inv(T);
    eCR = sol;
end