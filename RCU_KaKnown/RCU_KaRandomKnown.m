function [epsilon] = RCU_KaRandomKnown(P,P1,k,n,E_Ka,p_Ka) 
% function [epsilon] = RCU_KaRandomKnown(P,P1,k,n,E_Ka,p_Ka)
% Compute the RCU bound for the PUPE in Theorem 1 of
%
% [1] Y. Polyanskiy, "A perspective on massive random-access," 2017 IEEE 
% International Symposium on Information Theory (ISIT), 2017, pp.
% 2523-2527.
% 
% extended to the complex-valued model and the case where Ka is random and 
% unknown. See also
% 
% [2] K.-H. Ngo, A. Lancho, G. Durisi and A. G. i. Amat, "Unsourced
% Multiple Access With Random User Activity," submitted to IEEE Trans. Inf.
% Theory, Jan. 2022.
%
% INPUTS
%   P   : symbol power budget
%   P1  : the actual power used (denoted by P' in [1] and [2])
%   k   : number of bits per symbol
%   n   : framelength (number of complex DoFs)
%   E_Ka : average number of active users
%   p_Ka : PMF of the number of active users Ka
% 
% OUTPUTS
%   epsilon : (extended) bound on PUPE given in [1, Eq. (3)]

% codebook size
M       = 2^k;

% number of samples for empirical evaluation of the CDF of It in qt
Niq     = 1000;

%% Computation of p0 (with some additional terms w.r.t. p0 in [1] due to the randomness of Ka. See [2].)
% find K_l and K_u such that the probability that Ka is not in [K_l:K_u] is small 
K_l = floor(E_Ka); K_u = ceil(E_Ka);
while p_Ka(K_l-1) > 1e-9
    K_l = K_l - 1;
end
K_l = max(K_l,0);
while p_Ka(K_u+1) > 1e-9
    K_u = K_u + 1;
end

p01 = 0;
for Ka_tmp = K_l:K_u   
    p01_tmp = 1; 
    for ii = 1:Ka_tmp-1
        p01_tmp = p01_tmp*(1-ii/M);
    end
    p01 = p01 + p_Ka(Ka_tmp)*p01_tmp;
%     p01 = p01 + p_Ka(Ka_tmp)*nchoosek(Ka_tmp,2)/M;
end
p0 = 2 - p01 - sum(p_Ka(K_l:K_u)) + E_Ka*gammainc(n*P/P1,n,'upper');

%% Initialize epsilon to p0
epsilon = p0;

%% Compute the RCU bound, averaged over the distribution of Ka
parfor Ka = max(K_l,1):K_u  
    
    %% Computation of pt
    % Initialize for the sum over t and the optimization over rho and rho1
    rho_vec = linspace(0,1,100); 
    rho1_vec= linspace(0,1,100); 
    t_vec   = 1:Ka;
    rho     = repmat(rho_vec,length(rho1_vec),1,length(t_vec));
    rho1    = permute(repmat(rho1_vec,length(rho_vec),1,length(t_vec)),[2,1,3]);
    t       = permute(repmat(t_vec,length(rho1_vec),1,length(rho_vec)),[1 3 2]);

    % Some functions for the computation of pt as defined in [1, Th. 1]
    R1_f = @(n,M,t)  1./n*log(M)-1./(n*t).*gammaln(t+1);
    R2_f = @(n,Ka,t) 1/n*(gammaln(Ka+1)-gammaln(t+1)-gammaln(Ka-t+1));
    D_f  = @(P1,t,rho,rho1) (P1.*t-1).^2 + 4*P1.*t.*(1+rho.*rho1)./(1+rho);
    lambda_f = @(P1,t,rho,rho1) (P1.*t-1+sqrt(D_f(P1,t,rho,rho1)))./(2.*(1+rho1.*rho).*P1.*t); % max(roots_edit([P1.*t.*(rho+1).*(rho.*rho1+1) -P1.*t.*(rho+1) -1]));
    mu_f = @(P1,t,rho,rho1) rho.*lambda_f(P1,t,rho,rho1)./(1+P1.*t.*lambda_f(P1,t,rho,rho1));
    a_f = @(P1,t,rho,rho1) rho.*log(1+P1.*t.*lambda_f(P1,t,rho,rho1))+log(1+P1.*t.*mu_f(P1,t,rho,rho1));
    b_f = @(P1,t,rho,rho1) rho.*lambda_f(P1,t,rho,rho1)- mu_f(P1,t,rho,rho1)./(1+P1.*t.*mu_f(P1,t,rho,rho1));
    E0_f = @(P1,t,rho,rho1) rho1.*a_f(P1,t,rho,rho1) + log(1-b_f(P1,t,rho,rho1).*rho1);
    Et_f = @(P1,t,rho,rho1,n,M,Ka) squeeze(max(max(-rho.*rho1.*t.*R1_f(n,M,t) - rho1.*R2_f(n,Ka,t) + E0_f(P1,t,rho,rho1))));
    pt_f   = @(P1,t,rho,rho1,n,M,Ka) exp(-n*Et_f(P1,t,rho,rho1,n,M,Ka));

    % Compute pt
    pt = pt_f(P1,t,rho,rho1,n,M,Ka);

    %% Computation of qt (for t = 1 only, as in [1])
    qt = 1;
    % compute qt for Ka < 50 only due to complexity issue
    if Ka < 50 
        It = zeros(1,Niq);
        for II = 1:Niq
            Zi = sqrt(.5)*(randn(1,n) + 1i*randn(1,n));
            codebook = sqrt(.5*P1)*(randn(Ka,n) + 1i*randn(Ka,n));
            it   = n*log(1+P1) + (sum(abs(repmat(Zi,Ka,1)+codebook).^2,2)./(1+P1)-sum(abs(repmat(Zi,Ka,1)).^2,2));
            It(II) = min(it);
        end
        [prob,gam] = ecdf(It);
        qt = min(prob + exp(n*(R1_f(n,M,1)+R2_f(n,Ka,1))-gam));
    end
    
    %% Take the min of pt and qt 
    pt(1) = min(pt(1),qt);

    %% Compute epsilon, the RHS of [1, Eq. (3)]
    epsilon = epsilon + (t_vec/Ka*pt)*feval(p_Ka,Ka);
end
end