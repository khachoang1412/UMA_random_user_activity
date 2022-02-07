function [eps_MD,eps_FA] = RCU_KaRandomUnknown(P,P1,rad_l,rad_u,k,n,E_Ka,p_Ka)
% function [eps_MD,eps_FA] = RCU_KaRandomUnknown(P,P1,rad_l,rad_u,k,n,E_Ka,p_Ka)
% Compute the RCU bound for the PUPE in Theorem 1 of
%
% [1] K.-H. Ngo, A. Lancho, G. Durisi and A. G. i. Amat, "Unsourced
% Multiple Access With Random User Activity," submitted to IEEE Trans. Inf.
% Theory, Jan. 2022.
% 
% for a system with random and unknown number of active users.
%
% INPUTS
%   P   : symbol power budget
%   P1  : the actual power used (denoted by P' in [1] and [2])
%   rad_l : lower decoding radius
%   rad_u : upper decoding radius
%   k     : number of bits per symbol
%   n     : framelength (number of complex DoFs)
%   E_Ka  : average number of active users
%   p_Ka  : PMF of the number of active users Ka
% 
% OUTPUTS
%   eps_MD, eps_FA : RCU bounds on the misdetection and false alarm 
%                       probabilities, respectively

% codebook size
M       = 2^k;

% number of samples for empirical evaluation of the CDF of It in qt
Niq     = 1000;

%% Computation of \tilde{p}
% Compute the error floors [1, Th. 3]
[floor_MD,floor_FA] = RCU_floor_KaRandomUnknown(rad_l,rad_u,k,n,E_Ka,p_Ka);

% Find K_l and K_u such that Pr[Ka \notin [K_l:K_u]] is small
K_l = floor(E_Ka); K_u = ceil(E_Ka);
while p_Ka(K_l-1) > max(.0001*min(floor_MD,floor_FA),1e-9)
    K_l = K_l - 1;
end
K_l = max(K_l,0);
while p_Ka(K_u+1) > max(.0001*min(floor_MD,floor_FA),1e-9)
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
ptilde = 1 - p01 + 1 - sum(p_Ka(K_l:K_u)) + E_Ka*gammainc(n*P/P1,n,'upper');

if ptilde >= 1
    eps_MD = 1;
    eps_FA = 1;
    return
end

%% Initialize eps_MD and eps_FA to be p0
eps_MD = ptilde;
eps_FA = ptilde;

%% The expectation w.r.t. Ka
parfor Ka = K_l:K_u
%% Compute the term \xi(Ka,Ka') using [1, Th. 2]
% The function \zeta for ML esimator of Ka
KaEst_thres = @(Ka,KaEst,P1) n.*(log(1+Ka.*P1) - log(1+KaEst.*P1))./((1+Ka.*P1)./(1+KaEst.*P1)-1);

% For energy-based estimator of Ka, use the below
% KaEst_thres = @(Ka,KaEst,P1) n./(1+Ka.*P1).*(1+(Ka+KaEst).*P1/2); 

% Option 1: compute the exact term defined in [1]
% tic
% P_Ka_KaEst_list = zeros(Ka_u - Ka_l + 1,1);
% for idx = 1:Ka_u - Ka_l + 1
%     KaEst = Ka_l + idx - 1;
%     P_Ka_KaEst_list_a = gammainc(KaEst_thres(Ka,Ka_l:KaEst-1,P1),n,'lower');
%     P_Ka_KaEst_list_b = gammainc(KaEst_thres(Ka,KaEst+1:Ka_u,P1),n,'upper');
%     P_Ka_KaEst_list(idx) = 1 - max([P_Ka_KaEst_list_a P_Ka_KaEst_list_b]);
% end
% toc

% Option 2: slight relaxation that is faster
% tic
P_Ka_KaEst_list_a = gammainc(KaEst_thres(Ka,K_l:Ka-1,P1),n,'lower');
P_Ka_KaEst_list_b = gammainc(KaEst_thres(Ka,Ka+1:K_u,P1),n,'upper');
P_Ka_Ka = 1 - max([P_Ka_KaEst_list_a P_Ka_KaEst_list_b]);
P_Ka_KaEst_list = [P_Ka_KaEst_list_a'; P_Ka_Ka; P_Ka_KaEst_list_b'];
% toc

%% The expectation w.r.t. Ka'
for KaEst = K_l:K_u
    
    % \xi(Ka,Ka')
    xi_Ka_KaEst = P_Ka_KaEst_list(KaEst - K_l + 1);

    % \underbar{Ka'}
    KaEst_l = max(KaEst - rad_l,K_l);

    % \overbar{Ka'}
    KaEst_u = min(KaEst + rad_u,K_u);

    % Sets of possible values of t and t'
    t_vec = 0:min(min(KaEst_u,Ka),M-KaEst_l-max(Ka-KaEst_u,0));
    t1_l = max(Ka-KaEst_u,0) - max(Ka-KaEst_l,0);
    t1_u = max(KaEst_u-Ka,0) - max(KaEst_l-Ka,0);
    t1_vec = zeros(length(t_vec),KaEst_u-KaEst_l+1);
    for idxT = 1:length(t_vec)
        t1_vec(idxT,:) = t1_l+t_vec(idxT):t1_u+t_vec(idxT);
    end

    %% Computation of pt:

    % Definitions of functions R_1 and R_2
    if M-max(Ka,KaEst_l) <= 2^20
        R1_f = @(n,Ka,K1,M,t1) (gammaln(M-max(Ka,K1) + 1) - gammaln(M-max(Ka,K1)-t1+ 1) - gammaln(t1+1))/n;
    else
        R1_f = @(n,Ka,K1,M,t1) t1./n.*log(M-max(Ka,K1))-1./n.*gammaln(t1+1); 
    end
    R2_f = @(n,Ka,K1,t) 1/n*(gammaln(min(Ka,K1)+1)-gammaln(t+1)-gammaln(min(Ka,K1)-t+1));
    
    % First, consider t1 > 0
    t1_vec2 = t1_vec;
    t1_vec2(t1_vec2<0) = 0;
    
    % Initialize for the optimization over rho and rho1
    rho_vec =  linspace(1e-9,1,50); 
    rho1_vec = linspace(1e-9,1,50); 
    
    rho = permute(repmat(rho_vec,length(rho1_vec),1,length(t_vec),size(t1_vec,2)),[2,1,3,4]);
    rho1 = repmat(rho1_vec,length(rho_vec),1,length(t_vec),size(t1_vec,2));
    t = permute(repmat(t_vec,length(rho_vec),1,length(rho1_vec),size(t1_vec,2)),[1,3,2,4]);
    t1 = permute(repmat(t1_vec2,1,1,length(rho_vec),length(rho1_vec)),[3,4,1,2]);
    
    % Precompute some quantities
    P2 = 1+(max(Ka-KaEst_u,0)+max(KaEst_l-Ka,0)).*P1;
    P3 = P1.*(t1+rho.*t);
    P4 = rho.*rho1.*P3.*P2;
    P5 = rho1.*(rho-1).*t1.*P1;

    % Find optimal lambda that maximizes E0(t,t')
    c1_f = - P3.*P4.*P5 - (rho1+1).*t1.*P1.*P3.*P4;
    c2_f = P5.*(P3.^2 - P4) + rho1.*(t1.*P1.*P3.^2 - P3.*P4) ...
            - (2.*t1.*P1 + P3).*P4;
    c3_f = 2.*P3.*P5 + rho1.*(P3.^2 + t1.*P1.*P3) - 2.*P4;
    c4_f = P5 + rho1.*P3; 

    Delta0_f = c2_f.^2 - 3.*c1_f.*c3_f;
    Delta1_f = 2.*c2_f.^3 - 9.*c1_f.*c2_f.*c3_f + 27.*c1_f.^2.*c4_f;
    Q_f = ((Delta1_f + sqrt(Delta1_f.^2 - 4.*Delta0_f.^3))./2).^(1/3);

    lambda = real(-(c2_f + Q_f + Delta0_f./Q_f)./3./c1_f);
    
    % Compute E0(t,t') 
    E0_f = rho1.*(rho-1).*log(1+lambda.*P1.*t1) ...
        + (rho1-1).*log(1+lambda.*P3) ...
        + log(1+lambda.*P3 - lambda.^2.*P4);
    
    % Compute E(t,t') 
    Ett_f = reshape(max(max(-rho.*rho1.*R1_f(n,Ka,KaEst_l,M,t1) ...
        - rho1.*R2_f(n,Ka,KaEst_u,t) + E0_f)),[length(t_vec) size(t1_vec,2)]); 
    
    % Compute p_{t,t'}
    ptt = min(exp(-n.*Ett_f),1);

    % Now, consider the case t1 = 0, where the optimal lambda is the root 
    % of a quadratic function
    rho = permute(repmat(rho_vec,length(rho1_vec),1,length(t_vec)),[2,1,3]);
    rho1 = repmat(rho1_vec,length(rho_vec),1,length(t_vec));
    t = permute(repmat(t_vec,length(rho_vec),1,length(rho1_vec)),[1,3,2]);
    P3 = P1.*rho.*t;
    P4 = rho.*rho1.*P3.*P2;

    c2_f = -(rho1+1).*P3.*P4;
    c3_f = rho1.*P3.^2 - 2.*P4;
    c4_f = rho1.*P3;
    Delta = c3_f.^2 - 4.*c2_f.*c4_f;
    lambda = -(c3_f+sqrt(Delta))./2./c2_f;
    E0_f = (rho1-1).*log(1+lambda.*P3) + log(1+lambda.*P3 - lambda.^2.*P4);
    Ett_f = squeeze(max(max(- rho1.*R2_f(n,Ka,KaEst_u,t) + E0_f))); 
    ptt_t1zero = min(exp(-n.*Ett_f),1);
    
    % Combine the cases t1 > 0 and t1 = 0
    ptt(t1_vec > min(KaEst_u-max(KaEst_l-Ka,0),M-max(KaEst_l,Ka))) = 0;
    ptt(t1_vec < 0) = 0;
    for idxT = 1:length(t_vec)
    for idxT1 = 1:size(t1_vec,2)
        if t1_vec(idxT,idxT1) == 0
            ptt(idxT,idxT1) = ptt_t1zero(idxT);
        end
    end
    end

    % Compute pt
    pt = min(sum(ptt,2),1);
    
    %% Computation of qt for t = 1:
    qt = 1;
    qtt = 1;

    if t_vec(end) >= 1

    t1_set = t1_vec2(2,:);
    
    % Compute qt for t = 1 and Ka <= 50 only due to complexity issue
    if t_vec(1) <= 1 && Ka <= 150
        It = zeros(1,Niq);
        for II = 1:Niq
            Zi = sqrt(.5)*(randn(1,n) + 1i*randn(1,n));
            DefaultFA = sqrt(.5*max(KaEst_l-Ka,0)*P1)*(randn(1,n) + 1i*randn(1,n));
            codebook = sqrt(.5*P1)*(randn(min(Ka,KaEst_u),n) + 1i*randn(min(Ka,KaEst_u),n));

            it = n*log(1+(1+max(Ka-KaEst_u,0))*P1) + ...
                (sum(abs(repmat(Zi+DefaultFA,min(Ka,KaEst_u),1)+codebook).^2,2)./...
                    (1+(1+max(Ka-KaEst_u,0))*P1)-sum(abs(repmat(Zi,min(Ka,KaEst_u),1)).^2,2));
            It(II) = min(it);
        end
        [prob,gam] = ecdf(It);
        qt = min(prob + sum(exp(repmat(n*(R1_f(n,Ka,KaEst_l,M,t1_set)+R2_f(n,Ka,KaEst_u,1)),length(gam),1)-gam),2));
        qtt = min(repmat(prob,1,length(t1_set)) + exp(n*(R1_f(n,Ka,KaEst_l,M,t1_set)+R2_f(n,Ka,KaEst_u,1))-gam));
    end 
    
    % Take the min of pt and qt
    pt(2) = min(pt(2),qt);
    ptt(2,:) = min(ptt(2,:),qtt);
    end
    
    % Take the min of pt, qt, and xi(Ka,Ka')
    pt = min(pt,xi_Ka_KaEst);
    ptt = min(ptt,xi_Ka_KaEst);

    %% Computation epsilon_MD and epsilon_FA for a given P'< P
    if Ka > 0
        eps_MD = eps_MD + feval(p_Ka, Ka)*(sum((t_vec + max(Ka-KaEst_u,0))*pt)/Ka);
    end

    t_vec_2 = repmat(t_vec,size(t1_vec,2),1).';
    Mrx = (Ka-t_vec_2+t1_vec + max(KaEst_l-Ka,0) - max(Ka-KaEst_u,0)); % number of decoded codewords
    Mrx(Mrx == 0) = inf; 
    eps_FA = eps_FA + feval(p_Ka, Ka)*sum(sum(ptt.*(t1_vec + max(KaEst_l-Ka,0))./Mrx,1));
end
% if eps_MD >= 1 && eps_FA >= 1
%     break;
% end
end
eps_MD = min(eps_MD,1);
eps_FA = min(eps_FA,1);
end