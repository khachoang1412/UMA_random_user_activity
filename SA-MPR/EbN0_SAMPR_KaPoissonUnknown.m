function data = EbN0_SAMPR_KaPoissonUnknown(k, n, epsilon_MD, epsilon_FA, E_Ka, SlotIdxCoding)
% function data = EbN0_ALOHA_KaPoissonUnknown(k, n, epsilon, E_Ka, SlotIdxCoding)
% Find the minimal required EbN0 (in dB) such that the misdetection and 
% false alarm probabilities achieved with slotted ALOHA (SA) with 
% multi-packet reception (MPR) is below certain thresholds epsilon_MD and
% epsilon_FA, respectively, for a system with the number of active 
% users following a Poisson distribution, and unknown. See Corollary 2 in 
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
%   SlotIdxCoding : 1 if slot-index coding is employed, 0 otherwise 
%
% OUTPUT
%   data : store the system parameters and the minimal required EbN0

addpath ../RCU_KaUnknown

tic
DEBUG = 0;

%% debugging mode
if DEBUG == 1
    k       = 128; 
    n       = 19200; 
    epsilon_MD = .1; 
    epsilon_FA = .1; 
    E_Ka    = 50; 
    SlotIdxCoding = 0;
end

%% Poissom PMF of Ka, can be modified to consider other distributions 
p_Ka = @(K) poisspdf(K,E_Ka);

%% The range of power to search
EbN0db_lowest = -1;
EbN0db_highest = 16;
P_lowest = k*10^(EbN0db_lowest/10)/n;
P_highest = k*10^(EbN0db_highest/10)/n;

%% Function to compute the RCU bound
if SlotIdxCoding == 1
    f_rcu = @(P,P1,L,DecRad) RCU_KaRandomUnknown(P*L,P1*L,DecRad,DecRad,k-floor(log2(L)),...
                floor(n/L),E_Ka/L,@(Ksa) PMF_Ksa(p_Ka,E_Ka,L,Ksa));
else
    f_rcu = @(P,P1,L,DecRad) RCU_KaRandomUnknown(P*L,P1*L,DecRad,DecRad,k,floor(n/L),...
                E_Ka/L,@(Ksa) PMF_Ksa(p_Ka,E_Ka,L,Ksa));
end

%% Search for the minimal required EbN0
% The RCU is optimized over P1, L, and the decoding radius
[eps_RCU_MD, eps_RCU_FA, P_RCU,P1,L,DecRad] = binary_search_P_SAMPR(f_rcu, ...
    P_lowest, P_highest,epsilon_MD,epsilon_FA,min(epsilon_MD,epsilon_FA)/100,E_Ka);
EbN0db_RCU = 10*log10(n*P_RCU/k);

%% Save the results
sim_time = toc;
data.EbN0db = EbN0db_RCU;
data.E_Ka   = E_Ka;
data.p_Ka   = 'Poisson';
data.eps_RCU_MD = eps_RCU_MD;
data.eps_RCU_FA = eps_RCU_FA;
data.epsilon_MD = epsilon_MD;
data.epsilon_FA = epsilon_FA;
data.DecRad = DecRad;
data.k      = k;
data.n      = n;
data.P1     = P1;
data.nSlot  = L; 
data.SlotIdxCoding = SlotIdxCoding; 
data.sim_time = sim_time;

if DEBUG ~= 1
    if SlotIdxCoding == 1
        filename = ['EbN0_SAMPR_SlotIdxCoding_KaPoissonUnknown_EKa_' num2str(E_Ka)...
            '_epsilonMD_' num2str(epsilon_MD) '_epsilonFA_' num2str(epsilon_FA)...
                '_k_' num2str(k) '_n_' num2str(n) '.mat'];
    else
        filename = ['EbN0_SAMPR_KaPoissonUnknown_EKa_' num2str(E_Ka) ...
            '_epsilonMD_' num2str(epsilon_MD) '_epsilonFA_' num2str(epsilon_FA) ...
                '_k_' num2str(k) '_n_' num2str(n) '.mat'];
    end
    save(filename, 'data', '-v7.3');
else
    keyboard
end

end

%%
function PMF_Ksa = PMF_Ksa(p_Ka,E_Ka,L,Ksa)
% Compute the PMF of the number of active users per slot

PMF_Ksa = zeros(size(Ksa));
 
K_u = E_Ka;
while p_Ka(K_u+1) > eps
    K_u = K_u + 1;
end
K_u = max(K_u,10000);

PMF_Ksa_scalar = @(T) sum(p_Ka(0:K_u).*binopdf(T,0:K_u,1/L));
for idxKsa = 1:length(Ksa)
    PMF_Ksa(idxKsa) = PMF_Ksa_scalar(Ksa(idxKsa));
end
end

%%
function [rcu_MD, rcu_FA, P,P1,nSlot,DecRad] = binary_search_P_SAMPR(f,x1,x2,TARGET_MD,TARGET_FA,TOL,E_Ka)
% Search for the value of P such that the RCU bound of the PUPE of SAMPR 
% given by min_{P1,L,DecRad} weightedSum[f(P,P1,L,DecRad)] is between 
% TARGET - TOL and TARGET

weight_MD = 1/(1 + TARGET_MD/TARGET_FA);
weight_FA = 1/(1 + TARGET_FA/TARGET_MD);

iter= 20;                       % maximum number of iterations
k1=0;                            % number of iterations

[rcu_MD_tmp,rcu_FA_tmp,P1_tmp,nSlot_tmp,DecRad_tmp] = ...
    golden_search_P1_SAMPR(f,1e-9,x2,x2/100, weight_MD, weight_FA,E_Ka);
[rcu_MD_tmp1,rcu_FA_tmp1,P1_tmp1,nSlot_tmp1,DecRad_tmp1] = ...
    golden_search_P1_SAMPR(f,1e-9,x1,x1/100, weight_MD, weight_FA,E_Ka);
if x1 == x2
    P = x1;
    rcu_MD = rcu_MD_tmp;
    rcu_FA = rcu_FA_tmp;
    P1 = P1_tmp;
    nSlot = nSlot_tmp;
    DecRad = DecRad_tmp;
elseif rcu_MD_tmp >= TARGET_MD || rcu_FA_tmp >= TARGET_FA
    warning('Impossible to achieve the target within the given range of parameter :( ');
    P = x2;
    rcu_MD = rcu_MD_tmp;
    rcu_FA = rcu_FA_tmp;
    P1 = P1_tmp;
    nSlot = nSlot_tmp;
    DecRad = DecRad_tmp;
elseif rcu_MD_tmp1 <= TARGET_MD && rcu_FA_tmp1 <= TARGET_FA
    warning('All parameter values in the range can achieve the target :) ');
    P = x1;
    rcu_MD = rcu_MD_tmp1;
    rcu_FA = rcu_FA_tmp1;
    P1 = P1_tmp1;
    nSlot = nSlot_tmp1;
    DecRad = DecRad_tmp1;
else
    
[fx_MD,fx_FA,P1_tmp,nSlot_tmp,DecRad_tmp]=...
    golden_search_P1_SAMPR(f,1e-9,(x1+x2)/2,(x1+x2)/200, weight_MD, weight_FA,E_Ka);

while ~((TARGET_MD >= fx_MD && fx_MD >= TARGET_MD - TOL && TARGET_FA >= fx_FA) || ...
        (TARGET_FA >= fx_FA && fx_FA >= TARGET_FA - TOL && TARGET_MD >= fx_MD) ...
        || (k1>iter))
    if k1 > 0
        [fx_MD,fx_FA,P1_tmp,nSlot_tmp,DecRad_tmp]=...
            golden_search_P1_SAMPR(f,0,(x1+x2)/2,(x1+x2)/200, weight_MD, weight_FA,E_Ka); 
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
nSlot = nSlot_tmp;
DecRad = DecRad_tmp;
end
end

%%
function [rcu_MD, rcu_FA, P1,nSlot,DecRad] = golden_search_P1_SAMPR(f, START_INT, END_INT, TOL, weight_MD, weight_FA,E_Ka)
% Optimize P1

P = END_INT;
iter= 20;                       % maximum number of iterations
tau=double((sqrt(5)-1)/2);      % golden proportion coefficient, around 0.618
k2=0;                            % number of iterations
x1=START_INT+(1-tau)*(END_INT-START_INT);             % computing x values
x2=START_INT+tau*(END_INT-START_INT);

% K_l = Ka;
% K_u = 1.5*E_Ka;

TOL_nSlot = 1;
nSlot_l = E_Ka;
nSlot_u = 6*E_Ka;

[fMD_x1,fFA_x1,nSlot_1,DecRad_1] = ...
    golden_search_nSlot_SAMPR(f, P,x1,nSlot_l,nSlot_u,TOL_nSlot,weight_MD,weight_FA); 
[fMD_x2,fFA_x2,nSlot_2,DecRad_2] = ...
    golden_search_nSlot_SAMPR(f, P,x2,nSlot_l,nSlot_u,TOL_nSlot,weight_MD,weight_FA);

while ((abs(END_INT-START_INT)>TOL) && (k2<iter))
    if fMD_x1*weight_MD + fFA_x1*weight_FA < fMD_x2*weight_MD + fFA_x2*weight_FA %function smaller in x1 than in x2
        END_INT = x2; %make interval smaller
        x2 = x1; %set new end of interval
        x1 = START_INT+(1-tau)*(END_INT-START_INT); %find new beginning
        fMD_x2 = fMD_x1;%already have value in x1
        fFA_x2 = fFA_x1;
        DecRad_2 = DecRad_1;
        nSlot_2 = nSlot_1;
        [fMD_x1,fFA_x1,nSlot_1,DecRad_1] = ...
            golden_search_nSlot_SAMPR(f, P,x1,nSlot_l,nSlot_u,TOL_nSlot,weight_MD,weight_FA); %compute new value for new beginning
    else
        START_INT=x1; %function smaller in x2 than in x1 so set new start indx to x1
        x1 = x2; %replace as new start index
        x2=START_INT+tau*(END_INT-START_INT); %compute new end index
        fMD_x1= fMD_x2;
        fFA_x1 = fFA_x2;
        nSlot_1 = nSlot_2;
        DecRad_1 = DecRad_2;
        [fMD_x2,fFA_x2,nSlot_2,DecRad_2] = ...
            golden_search_nSlot_SAMPR(f, P,x2,nSlot_l,nSlot_u,TOL_nSlot,weight_MD,weight_FA); 
    end
    k2=k2+1;
end

if fMD_x1*weight_MD + fFA_x1*weight_FA < fMD_x2*weight_MD + fFA_x2*weight_FA 
    P1=x1;
    rcu_MD = fMD_x1;
    rcu_FA = fFA_x1;
    nSlot = nSlot_1;
    DecRad = DecRad_1;
else
    P1=x2;
    rcu_MD = fMD_x2;
    rcu_FA = fFA_x2;
    nSlot = nSlot_2;
    DecRad = DecRad_2;
end

end

%%
function [rcu_MD, rcu_FA, nSlot, DecRad] = golden_search_nSlot_SAMPR(f, P,P1,START_INT,END_INT,TOL,weight_MD,weight_FA)
% Optimize the number of slots

iter= 20;                       % maximum number of iterations
tau=double((sqrt(5)-1)/2);      % golden proportion coefficient, around 0.618
k2=0;                            % number of iterations
x1=max(round(START_INT+(1-tau)*(END_INT-START_INT)),0);             % computing x values
x2=round(START_INT+tau*(END_INT-START_INT));

[fMD_x1,fFA_x1,DecRad_x1] = linear_search_DecRad_SAMPR(f,P,P1,x1,weight_MD,weight_FA); % computing values in x points
[fMD_x2,fFA_x2,DecRad_x2] = linear_search_DecRad_SAMPR(f,P,P1,x2,weight_MD,weight_FA); 

while ((abs(END_INT-START_INT)>TOL) && (k2<iter))
    if fMD_x1*weight_MD + fFA_x1*weight_FA < fMD_x2*weight_MD + fFA_x2*weight_FA %function smaller in x1 than in x2
        END_INT = x2; %make interval smaller
        x2 = x1; %set new end of interval
        x1 = max(round(START_INT+(1-tau)*(END_INT-START_INT)),0); %find new beginning
        fMD_x2 = fMD_x1;%already have value in x1
        fFA_x2 = fFA_x1;
        DecRad_x2 = DecRad_x1;
        if fixListSize
            [fMD_x1,fFA_x1] = f(P,P1,x1); DecRad_x1 = 0;
        else
            [fMD_x1,fFA_x1,DecRad_x1] = linear_search_DecRad_SAMPR(f,P,P1,x1,weight_MD,weight_FA); %compute new value for new beginning
        end
    else
        START_INT=x1; %function smaller in x2 than in x1 so set new start indx to x1
        x1 = x2; %replace as new start index
        x2 = round(START_INT+tau*(END_INT-START_INT)); %compute new end index
        fMD_x1= fMD_x2;
        fFA_x1 = fFA_x2;
        DecRad_x2 = 0;
        if fixListSize
            [fMD_x2,fFA_x2] = f(P,P1,x2);
        else
            [fMD_x2,fFA_x2,DecRad_x2] = linear_search_DecRad_SAMPR(f,P,P1,x2,weight_MD,weight_FA); 
        end
    end
    k2=k2+1;
end

if fMD_x1*weight_MD + fFA_x1*weight_FA < fMD_x2*weight_MD + fFA_x2*weight_FA
    nSlot=x1;
    rcu_MD = fMD_x1;
    rcu_FA = fFA_x1;
    DecRad = DecRad_x1;
else
    nSlot=x2;
    rcu_MD = fMD_x2;
    rcu_FA = fFA_x2;
    DecRad = DecRad_x2;
end
end

%%
function [fMD,fFA,DecRad] = linear_search_DecRad_SAMPR(f,P,P1,x,weight_MD,weight_FA)
% Optimize the decoding radius

maxrad = 5; % search decoding radius from 0 to 4
fMD = zeros(maxrad,1);
fFA = zeros(maxrad,1);
for ii = 1:maxrad
    [fMD_tmp, fFA_tmp] = f(P,P1,x,ii-1);
    fMD(ii) = fMD_tmp;
    fFA(ii) = fFA_tmp;
end
[~,idxMin] = min(fMD*weight_MD + fFA*weight_FA);
fMD = fMD(idxMin);
fFA = fFA(idxMin);
iiset = 0:maxrad-1;
DecRad = iiset(idxMin);
% DecRad = 0;
% [fMD, fFA] = f(P,P1,x,DecRad);
end