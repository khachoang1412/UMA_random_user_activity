function data = RCU_KaPoissonUnknown(k, n, E_Ka_list, EbN0db, rad_l_list, rad_u_list)
% function data = RCU_KaPoissonUnknown(k, n, E_Ka_list, EbN0db, rad_l_list, rad_u_list)
% Compute the RCU bound for the PUPE for a system with number of users 
% unknown and following a Poission distribution. The bound is given in 
% Theorem 1 of
%
% [1] K.-H. Ngo, A. Lancho, G. Durisi and A. G. i. Amat, "Unsourced
% Multiple Access With Random User Activity," submitted to IEEE Trans. Inf.
% Theory, Jan. 2022.
%
% INPUTS
%   k   : number of bits per symbol
%   n   : framelength (number of complex DoFs)
%   E_Ka_list  : set of values for the average number of active users
%   EbN0db     : energy per bit in dB
%   rad_l_list : set of values for the lower decoding radius
%   rad_u_list : set of corresponding upper decoding radius
% Note: rad_l_list and rad_u_list are of the same length
%
% OUTPUTS
%   data : store the system parameters and the computed RCU bound

tic
DEBUG = 0;

%% debugging mode
if DEBUG == 1
    k       = 128; 
    n       = 19200; 
    EbN0db  = 4; 
    E_Ka_list   = 50;   
    rad_l_list  = [0];
    rad_u_list = [0];
end

%% codebook size
M = 2^k;

%% symbol power budget
P_list = k.*10.^(EbN0db./10)./n;

%% initialization
p_MD = ones(length(EbN0db),length(E_Ka_list),length(rad_l_list));
p_FA = ones(length(EbN0db),length(E_Ka_list),length(rad_l_list));
P1 = zeros(length(EbN0db),length(E_Ka_list),length(rad_l_list));

%% Compute the RCU
for idxEKa = 1:length(E_Ka_list)
E_Ka = E_Ka_list(idxEKa);
p_Ka    = @(K) poisspdf(K,E_Ka);
for idxRad = 1:length(rad_l_list)
rad_l = rad_l_list(idxRad); % lower decoding radius
rad_u = rad_u_list(idxRad); % upper decoding radius
for idxEbN0 = 1:length(EbN0db)
    P = P_list(idxEbN0);
    
    f_rcu = @(P,P1) RCU_KaRandomUnknown(P,P1,rad_l,rad_u,k,n,E_Ka,p_Ka);
    
    % Optimize P1 to minimize max(p_MD,p_FA)
    [p_MD_tmp,p_FA_tmp,P1_tmp] = golden_search_P1_MDFA(f_rcu,0,P,P/200,'max',[],[]); 
    
    % To minimize a weighted sum of p_MD and p_FA, use the below
%     weight_MD = 1;
%     weight_FA = 1;
%     [p_MD_tmp,p_FA_tmp,P1_tmp] = golden_search_P1_MDFA(f_rcu,0,P,P/100,'weighted',weight_MD,weight_FA); 
    
    p_MD(idxEbN0,idxEKa,idxRad) = p_MD_tmp;
    p_FA(idxEbN0,idxEKa,idxRad) = p_FA_tmp;
    P1(idxEbN0,idxEKa,idxRad) = P1_tmp;
end
end
end
p_MD = squeeze(p_MD);
p_FA = squeeze(p_FA);
P1 = squeeze(P1);

%% Save the results
sim_time = toc;
data.EbN0db = EbN0db;
data.E_Ka   = E_Ka_list;
data.p_Ka   = 'Poisson';
data.k      = k;
data.n      = n;
data.rad_lower = rad_l;
data.rad_upper = rad_u;
data.p_MD   = p_MD;
data.p_FA   = p_FA;
data.P1     = P1;
data.sim_time = sim_time;

if DEBUG ~= 1
        filename = ['RCU_KaPoissonUnknown_EKa_' num2str(min(E_Ka_list)) 'to' ...
            num2str(max(E_Ka_list)) '_k_' num2str(k) '_n_' num2str(n) ...
            '_radL_' sprintf('%d', rad_l_list) ...
            '_radU_' sprintf('%d', rad_u_list) ...
            '_EbN0db_' num2str(min(EbN0db)) 'to' num2str(max(EbN0db)) '.mat'];
    save(filename, 'data', '-v7.3');
else
    keyboard
end

end