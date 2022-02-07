function data = EbN0_KaPoissonUnknown(k, n, epsilon_MD, epsilon_FA, E_Ka, rad_l, rad_u)
% function data = EbN0_KaPoissonUnknown(k, n, epsilon_MD, epsilon_FA, E_Ka)
% Find the minimal required EbN0 (in dB) such that the misdetection and 
% false alarm probabilities are below the thresholds epsilon_MD and 
% epsilon_FA, respectively, for a system with the number of active 
% users following a Poisson distribution and unknown.
%
% INPUTS
%   k   : number of bits per symbol
%   n   : framelength (number of complex DoFs)
%   epsilon_MD : target misdetection probability
%   epsilon_FA : target false alarm probability
%   E_Ka  : average number of active users
%   rad_l : lower decoding radius
%   rad_u : upper decoding radius
%
% OUTPUT
%   data : store the system parameters and the minimal required EbN0

tic
DEBUG = 0;

%% debugging mode
if DEBUG == 1
    k       = 128; % Number of bits
    n       = 19200; % Frame length
    epsilon_MD = .1; % Per-user error probability
    epsilon_FA = .1; % Per-user error probability
    E_Ka    = 50;
    rad_l = 0;
    rad_u = 0;
end

%% Poissom PMF of Ka, can be modified to consider other distributions 
p_Ka    = @(K) poisspdf(K,E_Ka);

%% The range of power to search
EbN0db_lowest = -2;
EbN0db_highest = 15;
    
P_lowest = k*10^(EbN0db_lowest/10)/n;
P_highest = k*10^(EbN0db_highest/10)/n;

%% The below can be used to find suitable decoding radii
% rad_l = 0;  rad_u = 0;
% 
% [floor_MD,floor_FA] = RCU_floor(rad_l,rad_u,k,n,E_Ka,p_Ka);
% while floor_MD > epsilon_MD || floor_FA > epsilon_FA
%     if floor_MD > epsilon_MD
%         rad_u = rad_u + 1;
%     end
%     if floor_FA > epsilon_FA
%         rad_l = rad_l + 1;
%     end
%     [floor_MD,floor_FA] = RCU_floor(rad_l,rad_u,k,n,E_Ka,p_Ka);
% end
% if (epsilon_MD < 1e-1) || (epsilon_FA < 1e-1)
%     rad_u = rad_u + 2;
%     rad_l = rad_l + 2;
% end

%% Function to compute the RCU bound
f_rcu = @(P,P1) RCU_KaRandomUnknown(P,P1,rad_l,rad_u,k,n,E_Ka,p_Ka);

%% Search for the minimal required EbN0
[eps_RCU_MD, eps_RCU_FA, P_RCU,P1] = binary_search_P_MDFA(f_rcu, P_lowest, P_highest,...
    epsilon_MD,epsilon_FA,min(epsilon_MD,epsilon_FA)/100,0);
EbN0db_RCU = 10*log10(n*P_RCU/k);

%% Save the results
sim_time = toc;
data.EbN0db = EbN0db_RCU;
data.E_Ka   = E_Ka;
data.p_Ka   = 'Poisson';
data.eps_est_MD = eps_RCU_MD;
data.eps_est_FA = eps_RCU_FA;
data.epsilon_MD = epsilon_MD;
data.epsilon_FA = epsilon_FA;
data.rad_l = rad_l;
data.rad_u = rad_u;
data.k      = k;
data.n      = n;
data.P1     = P1;
data.sim_time = sim_time;

if DEBUG ~= 1
    filename = ['EbN0_KaPoissonUnknown_EKa_' num2str(E_Ka) '_epsilonMD_' ...
        num2str(epsilon_MD) '_epsilonFA_' num2str(epsilon_FA) ...
        '_k_' num2str(k) '_n_' num2str(n) '.mat'];
    save(filename, 'data', '-v7.3');
else
    keyboard
end

end