function data = EbN0_KaPoissonKnown(k, n, epsilon, E_Ka)
% function data = EbN0_KaPoissonKnown(k, n, epsilon, E_Ka)
% Find the minimal required EbN0 (in dB) such that the PUPE is below a 
% certain threshold epsilon for a system with the number of active 
% users following a Poisson distribution but known.
%
% INPUTS
%   k   : number of bits per symbol
%   n   : framelength (number of complex DoFs)
%   epsilon : target PUPE
%   E_Ka    : average number of active users
%
% OUTPUT
%   data : store the system parameters and the minimal required EbN0

tic
DEBUG = 1;

%% debugging mode
if DEBUG == 1 
    k       = 100; 
    n       = 15000;
    epsilon = 0.1; 
    E_Ka      = 2; 
end

%% Poissom PMF of Ka, can be modified to consider other distributions 
p_Ka = @(K) poisspdf(K,E_Ka); 

%% The range of power to search
EbN0db_lowest = -0.5;
EbN0db_highest = 3;
    
P_lowest = k*10^(EbN0db_lowest/10)/n;
P_highest = k*10^(EbN0db_highest/10)/n;

%% Function to compute the RCU bound
f_rcu = @(P,P1) RCU_KaRandomKnown(P,P1,k,n,E_Ka,p_Ka);

%% Search for the minimal required EbN0
[eps_RCU, P_RCU] = binary_search(f_rcu, P_lowest, P_highest, epsilon,epsilon/100);
EbN0db_RCU = 10*log10(n*P_RCU/(k));

%% Save the results
sim_time = toc;
data.EbN0db = EbN0db_RCU;
data.E_Ka   = E_Ka;
data.eps_est= eps_RCU;
data.epsilon= epsilon;
data.k      = k;
data.n      = n;
data.sim_time = sim_time;

if DEBUG ~= 1
    filename = ['EbN0_KaPoissonKnown_EKa_' num2str(E_Ka) '_epsilon_' num2str(epsilon) '_k_' num2str(k) '_n_' num2str(n) '.mat'];
    save(filename, 'data', '-v7.3');
else
    keyboard
end

end