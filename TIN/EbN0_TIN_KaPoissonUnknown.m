function data = EbN0_TIN_KaPoissonUnknown(k, n, epsilon_MD, epsilon_FA, E_Ka, normalApprox)
% function data = EbN0_TIN_KaPoissonUnknown(k, n, epsilon_MD, epsilon_FA, E_Ka, normalApprox)
% Find the minimal required EbN0 (in dB) such that the misdetection and 
% false alarm probabilities achieved with treating-interference-as-noise
% (TIN) is below certain thresholds epsilon_MD and
% epsilon_FA, respectively, for a system with the number of active 
% users following a Poisson distribution, and unknown. See Theorem 4 in 
%
% [1] K.-H. Ngo, A. Lancho, G. Durisi and A. G. i. Amat, "Unsourced
% Multiple Access With Random User Activity," submitted to IEEE Trans. Inf.
% Theory, Jan. 2022.
%
% INPUTS
%   k   : number of bits per symbol
%   n   : framelength (number of complex DoFs)
%   epsilon_MD : target misdetection probability
%   epsilon_FA : target false alarm probability
%   E_Ka    : average number of active users
%   normalApprox : 1 if normal approximation is used, 0 otherwise 
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
    E_Ka    = 100; 
    normalApprox = 1;
end

%% Poissom PMF of Ka, can be modified to consider other distributions 
p_Ka    = @(K) poisspdf(K,E_Ka);

%% The range of power to search
EbN0db_lowest = 0;
EbN0db_highest = 15;
    
P_lowest = k*10^(EbN0db_lowest/10)/n;
P_highest = k*10^(EbN0db_highest/10)/n;

%% Function to compute the RCUs-TIN bound
f_TIN = @(P,P1,s) RCUs_TIN(P,P1,s,k,n,E_Ka,p_Ka,normalApprox);

%% Search for the minimal required EbN0
[eps_TIN_MD, eps_TIN_FA, P,P1,sopt] = binary_search_TIN(f_TIN, P_lowest, P_highest,...
    epsilon_MD,epsilon_FA,min(epsilon_MD,epsilon_FA)/100,normalApprox);
EbN0db = 10*log10(n*P/k);

%% Save the results
sim_time = toc;
data.EbN0db = EbN0db;
data.E_Ka   = E_Ka;
data.p_Ka   = 'Poisson';
data.eps_TIN_MD = eps_TIN_MD;
data.eps_TIN_FA = eps_TIN_FA;
data.epsilon_MD = epsilon_MD;
data.epsilon_FA = epsilon_FA;
data.k      = k;
data.n      = n;
data.P1     = P1;
data.s      = sopt;
data.normalApprox = normalApprox;
data.sim_time = sim_time;

if DEBUG ~= 1
    filename = ['EbN0_TIN_KaPoissonUnknown_EKa_' num2str(E_Ka) '_epsilonMD_' ...
        num2str(epsilon_MD) '_epsilonFA_' num2str(epsilon_FA) ...
        '_k_' num2str(k) '_n_' num2str(n) '.mat'];
    save(filename, 'data', '-v7.3');
else
    keyboard
end

end

%%
function Pe = eta_s(n,P,R,Ka_vec,s,N)
% Computation \eta_s in [1, th. 4], which is also a RCUs bound for the 
% noncoherent Rayleigh block-fading channel.
%
% INPUTS:
% n     = blocklength
% P     = transmit power (no dB!)
% R     = Values of Rate ln(M)/n at which the bound will be computed
% Ka    = Number of active users
% s     = parameter to optimize in the RCUs (s>0)
% N     = Number of samples to generate the generalized info. density
%
% OUTPUT:
% Pe = Error probability obtained as a result of the computation of the bound

Pe  = nan(size(Ka_vec));
parfor ii = 1:length(Ka_vec)
    Ka = Ka_vec(ii);
    snr = P/(1+(Ka-1)*P);
    x   = sqrt(0.5*snr)*(randn(n,N) + 1i*randn(n,N)); % signal tx by one of the users
    z   = sqrt(0.5)*(randn(n,N) + 1i*randn(n,N)); % Noise plus interference
    y   = x + z; % Received signal
    
    i_s = -s*vecnorm(y-x).^2 + s*vecnorm(y).^2/(1+s*snr) + n*log(1+s*snr); % generalized information density
    
    Pe(ii)  = mean( exp( - max(0, i_s - log(exp(n*R)-1))));
end
end

%%
function [eps_MD,eps_FA] = RCUs_TIN(P,P1,s,k,n,E_Ka,p_Ka,normalApprox)

% codebook size
M       = 2^k;

%% COMPUTATION OF \tilde{p}
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
end
ptilde = E_Ka*gammainc(n*P/P1,n,'upper') + 1 - p01 + 1 - sum(p_Ka(K_l:K_u)); 

%% Initialize eps_MD and eps_FA to be \tilde{p}
eps_MD = ptilde;
eps_FA = ptilde;

%% Compute the \xi(Ka,Ka') term for Ka and Ka' from K_l to K_u
if K_l == K_u
    xi_Ka_KaEst = 1;
else
    xi_Ka_KaEst = zeros(K_u-K_l+1,K_u-K_l+1);
    for Ka = K_l:K_u
        KaEst_thres = @(Ka,KaEst,P1) n.*(log(1+Ka.*P1) - log(1+KaEst.*P1))./((1+Ka.*P1)./(1+KaEst.*P1)-1);
        P_Ka_KaEst_list_a = gammainc(KaEst_thres(Ka,K_l:Ka-1,P1),n,'lower');
        P_Ka_KaEst_list_b = gammainc(KaEst_thres(Ka,Ka+1:K_u,P1),n,'upper');
        P_Ka_Ka = 1 - max([P_Ka_KaEst_list_a P_Ka_KaEst_list_b]);
        xi_Ka_KaEst(:,Ka-K_l+1) = [P_Ka_KaEst_list_a'; P_Ka_Ka; P_Ka_KaEst_list_b'];
    end
end

%% Compute the \eta_s terms
N = 1e2;
R = k/n*log(2);

if normalApprox
    C_TIN = @(Ka,P) log(1+P./(1+P.*(Ka-1)));
    V_TIN = @(Ka,P) sqrt(2.*P.*(log(exp(1))).^2./(1+Ka.*P));
    eta_s_vec = qfunc((n.*C_TIN((K_l:K_u)',P1) - log(M-1)) ./ (sqrt(n).*V_TIN((K_l:K_u)',P1)));
else
    eta_s_vec = eta_s(n,P1,R,(K_l:K_u)',s,N);
end

%% The sum over Ka
parfor Ka = K_l:K_u
    nAdditionalMDFA = min(min((K_l:K_u)',M-max(Ka,(K_l:K_u)')),Ka);
    eta = min(eta_s_vec,xi_Ka_KaEst(Ka-K_l+1,:)');
    
    % the sum over Ka' is implemented by the inner products 
    if Ka > 0
        eps_MD = eps_MD + (feval(p_Ka,Ka)/Ka)*(nAdditionalMDFA'*eta + max(Ka-(K_l:K_u),0)*xi_Ka_KaEst(Ka-K_l+1,:)');
    end
    if K_l == 0
        eps_FA = eps_FA + ...
            feval(p_Ka,Ka)*sum((nAdditionalMDFA(2:end).*eta(2:end) + max((K_l+1:K_u)'-Ka,0).*xi_Ka_KaEst(Ka-K_l+1,2:end)')./(K_l+1:K_u)');
    else
        eps_FA = eps_FA + ...
            feval(p_Ka,Ka)*sum((nAdditionalMDFA.*eta + max((K_l:K_u)'-Ka,0).*xi_Ka_KaEst(Ka-K_l+1,:)')./(K_l:K_u)');
    end
    
%     if eps_MD >= 1 && eps_FA >= 1
%         break;
%     end
end
eps_MD = min(eps_MD,1);
eps_FA = min(eps_FA,1);
end

%%
function [rcu_MD, rcu_FA, P,P1,s_opt] = binary_search_TIN(f,x1,x2,TARGET_MD,TARGET_FA,TOL,normalApprox)
% Search for P

weight_MD = 1/(1 + TARGET_MD/TARGET_FA);
weight_FA = 1/(1 + TARGET_FA/TARGET_MD);

iter= 20;                       % maximum number of iterations
k1=0;                            % number of iterations

[rcu_MD_tmp,rcu_FA_tmp,P1_tmp,s_tmp] = golden_search_P1_TIN(f,1e-9,x2,x2/200, weight_MD, weight_FA,normalApprox);
[rcu_MD_tmp1,rcu_FA_tmp1,P1_tmp1,s_tmp1] = golden_search_P1_TIN(f,1e-9,x1,x1/200, weight_MD, weight_FA,normalApprox);
if x1 == x2
    P = x1;
    rcu_MD = rcu_MD_tmp;
    rcu_FA = rcu_FA_tmp;
    P1 = P1_tmp;
    s_opt = s_tmp;
elseif rcu_MD_tmp > TARGET_MD || rcu_FA_tmp > TARGET_FA
    warning('Impossible to achieve the target within the given range of parameter :( ');
    P = x2;
    rcu_MD = rcu_MD_tmp;
    rcu_FA = rcu_FA_tmp;
    P1 = P1_tmp;
    s_opt = s_tmp1;
elseif rcu_MD_tmp1 < TARGET_MD || rcu_FA_tmp1 < TARGET_FA
    warning('All parameter values in the range can achieve the target :) ');
    P = x1;
    rcu_MD = rcu_MD_tmp1;
    rcu_FA = rcu_FA_tmp1;
    P1 = P1_tmp1;
    s_opt = s_tmp1;
else
    
[fx_MD,fx_FA,P1_tmp,s_tmp]=golden_search_P1_TIN(f,1e-9,(x1+x2)/2,(x1+x2)/200, weight_MD, weight_FA,normalApprox); % computing values in x points

while ~((TARGET_MD >= fx_MD && fx_MD >= TARGET_MD - TOL && TARGET_FA >= fx_FA) || ...
        (TARGET_FA >= fx_FA && fx_FA >= TARGET_FA - TOL && TARGET_MD >= fx_MD) ...
        || (k1>iter)) 
    if k1 > 0
        [fx_MD,fx_FA,P1_tmp,s_tmp]=golden_search_P1_TIN(f,0,(x1+x2)/2,(x1+x2)/200, weight_MD, weight_FA,normalApprox); 
    end
    if TARGET_MD > fx_MD && TARGET_FA > fx_FA
        x2 = (x1+x2)/2; %set new end of interval        
    else
        x1 = (x1+x2)/2; %replace as new start index
    end
    k1=k1+1;
end

rcu_MD = fx_MD;
rcu_FA = fx_FA;
P   = x2;
P1 = P1_tmp;
s_opt = s_tmp;
end
end

%%
function [rcu_MD, rcu_FA, P1,s_opt] = golden_search_P1_TIN(f, START_INT, END_INT, TOL, weight_MD, weight_FA,normalApprox)
% Optimize P1

P = END_INT;
iter= 20;                       % maximum number of iterations
tau=double((sqrt(5)-1)/2);      % golden proportion coefficient, around 0.618
k2=0;                            % number of iterations
x1=START_INT+(1-tau)*(END_INT-START_INT);             % computing x values
x2=START_INT+tau*(END_INT-START_INT);

s_start = 0;
s_end = 2;

if normalApprox
    [fMD_x1,fFA_x1] = f(P,x1,1);    s_opt1 = 1;
    [fMD_x2,fFA_x2] = f(P,x2,1);    s_opt2 = 1;
else
    [fMD_x1,fFA_x1,s_opt1] = golden_search_s_TIN(f, P, x1, s_start, s_end, TOL, weight_MD, weight_FA);
    [fMD_x2,fFA_x2,s_opt2] = golden_search_s_TIN(f, P, x2, s_start, s_end, TOL, weight_MD, weight_FA);
end

while ((abs(END_INT-START_INT)>TOL) && (k2<iter))
    if fMD_x1*weight_MD + fFA_x1*weight_FA < fMD_x2*weight_MD + fFA_x2*weight_FA %function smaller in x1 than in x2
        END_INT = x2; %make interval smaller
        x2 = x1; %set new end of interval
        x1 = START_INT+(1-tau)*(END_INT-START_INT); %find new beginning
        fMD_x2 = fMD_x1; %already have value in x1
        fFA_x2 = fFA_x1;
        s_opt2 = s_opt1;
        if normalApprox
            [fMD_x1,fFA_x1] = f(P,x1,1);    s_opt1 = 1;
        else
            [fMD_x1,fFA_x1,s_opt1] = golden_search_s_TIN(f, P, x1, s_start, s_end, TOL, weight_MD, weight_FA); %%compute new value for new beginning
        end
    else
        START_INT=x1; %function smaller in x2 than in x1 so set new start indx to x1
        x1 = x2; %replace as new start index
        x2=START_INT+tau*(END_INT-START_INT); %compute new end index
        fMD_x1= fMD_x2;
        fFA_x1 = fFA_x2;
        s_opt1 = s_opt2;
        if normalApprox
            [fMD_x2,fFA_x2] = f(P,x2,1);    s_opt2 = 1;
        else
            [fMD_x2,fFA_x2,s_opt2] = golden_search_s_TIN(f, P, x2, s_start, s_end, TOL, weight_MD, weight_FA);
        end
    end
    k2=k2+1;
end

if fMD_x1*weight_MD + fFA_x1*weight_FA < fMD_x2*weight_MD + fFA_x2*weight_FA 
    P1=x1;
    s_opt = s_opt1;
    rcu_MD = fMD_x1;
    rcu_FA = fFA_x1;
else
    P1=x2;
    s_opt = s_opt2;
    rcu_MD = fMD_x2;
    rcu_FA = fFA_x2;
end

end

%%
function [rcu_MD, rcu_FA, s_opt] = golden_search_s_TIN(f, P, P1, START_INT, END_INT, TOL, weight_MD, weight_FA)
% Optimize s

iter= 20;                       % maximum number of iterations
tau=double((sqrt(5)-1)/2);      % golden proportion coefficient, around 0.618
k2=0;                            % number of iterations
x1=START_INT+(1-tau)*(END_INT-START_INT);             % computing x values
x2=START_INT+tau*(END_INT-START_INT);

[fMD_x1,fFA_x1] = f(P,P1,x1);
[fMD_x2,fFA_x2] = f(P,P1,x2);

while ((abs(END_INT-START_INT)>TOL) && (k2<iter))
    if fMD_x1*weight_MD + fFA_x1*weight_FA < fMD_x2*weight_MD + fFA_x2*weight_FA %function smaller in x1 than in x2
        END_INT = x2; %make interval smaller
        x2 = x1; %set new end of interval
        x1 = START_INT+(1-tau)*(END_INT-START_INT); %find new beginning
        fMD_x2 = fMD_x1;%already have value in x1
        fFA_x2 = fFA_x1;
        [fMD_x1,fFA_x1] = f(P,P1,x1); %%compute new value for new beginning
    else
        START_INT=x1; %function smaller in x2 than in x1 so set new start indx to x1
        x1 = x2; %replace as new start index
        x2=START_INT+tau*(END_INT-START_INT); %compute new end index
        fMD_x1= fMD_x2;
        fFA_x1 = fFA_x2;
        [fMD_x2,fFA_x2] = f(P,P1,x2);
    end
    k2=k2+1;
end

if fMD_x1*weight_MD + fFA_x1*weight_FA < fMD_x2*weight_MD + fFA_x2*weight_FA 
    s_opt=x1;
    rcu_MD = fMD_x1;
    rcu_FA = fFA_x1;
else
    s_opt=x2;
    rcu_MD = fMD_x2;
    rcu_FA = fFA_x2;
end

end