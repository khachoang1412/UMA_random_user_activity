function [floor_MD,floor_FA] = RCU_floor_KaRandomUnknown(rad_l,rad_u,k,n,E_Ka,p_Ka)
% function [floor_MD,floor_FA] = RCU_floor(rad_l,rad_u,k,n,E_Ka,p_Ka)
% Compute the error floors for the RCU bounds of the misdetection and false
% alarm probabilities, characterized in Theorem 3 of
%
% [1] K.-H. Ngo, A. Lancho, G. Durisi and A. G. i. Amat, "Unsourced
% Multiple Access With Random User Activity," submitted to IEEE Trans. Inf.
% Theory, Jan. 2022.
% 
% for a system with random and unknown number of active users.
%
% INPUTS
%   rad_l : lower decoding radius
%   rad_u : upper decoding radius
%   k     : number of bits per symbol
%   n     : framelength (number of complex DoFs)
%   E_Ka  : average number of active users
%   p_Ka  : PMF of the number of active users Ka
% 
% OUTPUTS
%   floor_MD, floor_FA : the error floors

% codebook size
M       = 2^k;

%% Computation of \bar{p}
K_l = floor(E_Ka); K_u = ceil(E_Ka);
while p_Ka(K_l-1) >  1e-12
    K_l = K_l - 1;
end
K_l = max(K_l,0);
while p_Ka(K_u+1) > 1e-12
    K_u = K_u + 1;
end

p01 = 0;
for Ka_tmp = K_l:K_u   
    p01_tmp = 1; 
    for ii = 1:Ka_tmp-1
        p01_tmp = p01_tmp*(1-ii/M);
    end
    p01 = p01 + p_Ka(Ka_tmp)*p01_tmp;
end
pbar = 1 - p01 + 1 - sum(p_Ka(K_l:K_u));

%% Initialize floor_MD and floor_FA to be \bar{p}
floor_MD = pbar;
floor_FA = pbar;

%% The sum over Ka
parfor Ka = K_l:K_u  
% Compute \xi(Ka,Ka')
KaEst_thres = @(Ka,KaEst) n.*log(Ka./KaEst)./(Ka./KaEst-1);

P_Ka_KaEst_list_a = gammainc(KaEst_thres(Ka,K_l:Ka-1),n,'lower');
P_Ka_KaEst_list_b = gammainc(KaEst_thres(Ka,Ka+1:K_u),n,'upper');
P_Ka_Ka = 1 - max([P_Ka_KaEst_list_a P_Ka_KaEst_list_b]);
P_Ka_KaEst_list = [P_Ka_KaEst_list_a'; P_Ka_Ka; P_Ka_KaEst_list_b'];

%% The sum over Ka'
for KaEst = K_l:K_u
    if KaEst == 0
        P_Ka_KaEst = 0;
    else
        P_Ka_KaEst = P_Ka_KaEst_list(KaEst - K_l + 1);
    end

    KaEst_l = max(KaEst - rad_l,K_l);
    KaEst_u = min(KaEst + rad_u,K_u);

    if Ka > 0
        floor_MD = floor_MD + feval(p_Ka, Ka)*max(Ka-KaEst_u,0)*P_Ka_KaEst/Ka;
    end

    Mrx = (Ka + max(KaEst_l-Ka,0) - max(Ka-KaEst_u,0)); % number of decoded codewords
    if Mrx > 0
        floor_FA = floor_FA + feval(p_Ka, Ka)*max(KaEst_l-Ka,0)*P_Ka_KaEst/Mrx;
    end
end
end
floor_MD = min(floor_MD,1);
floor_FA = min(floor_FA,1);
end