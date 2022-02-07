function data = RCU_KaPoissonKnown(k, n, E_Ka, EbN0db)
% function data = RCU_KaPoissonKnown(k, n, E_Ka, EbN0db)
% Compute the RCU bound for the PUPE in Theorem 1 of
%
% [1] Y. Polyanskiy, "A perspective on massive random-access," 2017 IEEE 
% International Symposium on Information Theory (ISIT), 2017, pp.
% 2523-2527.
% 
% extended to the complex-valued model and the case where Ka follows a 
% Poisson distribution and is unknown. See also
% 
% [2] K.-H. Ngo, A. Lancho, G. Durisi and A. G. i. Amat, "Unsourced
% Multiple Access With Random User Activity," submitted to IEEE Trans. Inf.
% Theory, Jan. 2022.
%
% INPUTS
%   k   : number of bits per symbol
%   n   : framelength (number of complex DoFs)
%   E_Ka   : average number of active users
%   EbN0db : energy per bit in dB
% Note: E_Ka and EbN0db can be vectors.
%
% OUTPUTS
%   data : store the system parameters and the computed RCU bound

%   P   : symbol power budget
%   P1  : the actual power used (denoted by P' in [1] and [2])

tic
DEBUG = 0;

%% debugging mode
if DEBUG == 1
    k    = 128; 
    n    = 19200; 
    E_Ka = 50; 
    EbN0db = 0; 
end

%% symbol power budget
P_list = k.*10.^(EbN0db./10)./n;

%% initialization
epsilon = zeros(length(EbN0db),length(E_Ka));
P1 = zeros(length(EbN0db),length(E_Ka));

%% Compute the RCU bound
for idxEKa = 1:length(E_Ka)
    E_Ka = E_Ka(idxEKa);
    p_Ka = @(K) poisspdf(K,E_Ka);
    for idxEbN0 = 1:length(EbN0db)
        P = P_list(idxEbN0);

        f_rcu = @(P,P1) RCU(P,P1,k,n,E_Ka,p_Ka);
        
        [Pe_tmp,P1_tmp] = golden_search(f_rcu,0,P,P/100); % Optimize P'

        epsilon(idxEbN0,idxEKa) = Pe_tmp;
        P1(idxEbN0,idxEKa) = P1_tmp;
    end
end

%% Save the results
sim_time = toc;
data.EbN0db = EbN0db;
data.E_Ka    = E_Ka;
data.p_Ka   = 'Poisson';
data.epsilon= epsilon;
data.k      = k;
data.n      = n;
data.P1     = P1;
data.sim_time = sim_time;

if DEBUG ~= 1
    filename = ['RCU_KaPoissonKnown_EKa_' num2str(min(E_Ka)) 'to' ...
            num2str(max(E_Ka)) '_k_' num2str(k) '_n_' num2str(n) ...
            '_EbN0db_' num2str(min(EbN0db)) 'to' num2str(max(EbN0db)) '.mat'];
    save(filename, 'data', '-v7.3');
else
    keyboard
end

end