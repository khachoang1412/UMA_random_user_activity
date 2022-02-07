function data = EbN0_SAMPR_KaPoissonKnown(k, n, epsilon, E_Ka, SlotIdxCoding)
% function data = EbN0_SAMPR_KaPoissonKnown(k, n, epsilon, E_Ka, SlotIdxCoding)
% Find the minimal required EbN0 (in dB) such that the PUPE achieved with 
% slotted ALOHA (SA) with multi-packet reception (MPR) is below a 
% certain threshold epsilon for a system with the number of active 
% users following a Poisson distribution, but known.
%
% INPUTS
%   k   : number of bits per symbol
%   n   : framelength (number of complex DoFs)
%   epsilon : target PUPE
%   E_Ka    : average number of active users
%   SlotIdxCoding : 1 if slot-index coding is employed, 0 otherwise
%
% OUTPUT
%   data : store the system parameters and the minimal required EbN0

addpath ../RCU_KaKnown

tic
DEBUG = 0;

%% debugging mode
if DEBUG == 1
    k       = 128; 
    n       = 19200; 
    epsilon = .1; 
    E_Ka    = 50;
end

%% Poissom PMF of Ka, can be modified to consider other distributions 
p_Ka    = @(K) poisspdf(K,E_Ka);

%% The range of power to search
EbN0db_lowest = -1;
EbN0db_highest = 15;
    
P_lowest = k*10^(EbN0db_lowest/10)/n;
P_highest = k*10^(EbN0db_highest/10)/n;

%% Function to compute the RCU bound
if SlotIdxCoding == 1
    f_rcu = @(P,P1,L) RCU_KaRandomKnown(P*L,P1*L,k-floor(log2(L)),...
                ceil(n/L),E_Ka/L,@(Ksa) PMF_Ksa(p_Ka,E_Ka,L,Ksa)); 
else
    f_rcu = @(P,P1,L) RCU_KaRandomKnown(P*L,P1*L,k,ceil(n/L),E_Ka/L,...
                @(Ksa) PMF_Ksa(p_Ka,E_Ka,L,Ksa)); 
end

%% Search for the minimal required EbN0
% The RCU is optimized over P1 and L
[eps_RCU,P_RCU,P1,L] = binary_search_P_SAMPR_KaKnown(f_rcu, P_lowest, P_highest,...
    epsilon,epsilon/100,E_Ka);
EbN0db_RCU = 10*log10(n*P_RCU/k);

%% Save the results
sim_time = toc;
data.EbN0db = EbN0db_RCU;
data.E_Ka   = E_Ka;
data.p_Ka   = 'Poisson';
data.eps_RCU = eps_RCU;
data.epsilon = epsilon;
data.k      = k;
data.n      = n;
data.P1     = P1;
data.nSlot  = L; 
data.SlotIdxCoding  = SlotIdxCoding; 
data.sim_time = sim_time;

if DEBUG ~= 1
    if SlotIdxCoding == 1
        filename = ['EbN0_SAMPR_SlotIdxCoding_KaPoissonKnown_EKa_' num2str(E_Ka) ...
            '_epsilon_' num2str(epsilon) '_k_' num2str(k) '_n_' num2str(n) '.mat'];
    else
        filename = ['EbN0_SAMPR_KaPoissonKnown_EKa_' num2str(E_Ka) ...
            '_epsilon_' num2str(epsilon) '_k_' num2str(k) '_n_' num2str(n) '.mat'];
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
function [rcu,P,P1,nSlot] = binary_search_P_SAMPR_KaKnown(f,x1,x2,TARGET,TOL,E_Ka)
% Search for the value of P such that the RCU bound of the PUPE of SAMPR 
% given by min_{P1,L} f(P,P1,L) is between TARGET - TOL and TARGET

iter= 20;                       % maximum number of iterations
k1=0;                            % iteration index

[rcu_tmp,P1_tmp,nSlot_tmp] = golden_search_P1_SAMPR_KaKnown(f,1e-9,x2,x2/100,E_Ka);
[rcu_tmp1,P1_tmp1,nSlot_tmp1] = golden_search_P1_SAMPR_KaKnown(f,1e-9,x1,x1/100,E_Ka);

if x1 == x2
    P = x1;
    rcu = rcu_tmp;
    P1 = P1_tmp;
    nSlot = nSlot_tmp;
elseif rcu_tmp >= TARGET
    warning('Impossible to achieve the target within the given range of parameter :( ');
    P = x2;
    rcu = rcu_tmp;
    P1 = P1_tmp;
    nSlot = nSlot_tmp;
elseif rcu_tmp1 <= TARGET
    warning('All parameter values in the range can achieve the target :) ');
    P = x1;
    rcu = rcu_tmp1;
    P1 = P1_tmp1;
    nSlot = nSlot_tmp1;
else
    
[fx,P1_tmp,nSlot_tmp]=...
    golden_search_P1_SAMPR_KaKnown(f,1e-9,(x1+x2)/2,(x1+x2)/200,E_Ka); 

while (TARGET < fx || fx < TARGET - TOL) && (k1<iter || TARGET < fx)
    if k1 > 0
        [fx,P1_tmp,nSlot_tmp]=...
            golden_search_P1_SAMPR_KaKnown(f,0,(x1+x2)/2,(x1+x2)/200,E_Ka); 
    end
    if TARGET > fx
        x2 = (x1+x2)/2; %set new end of interval        
    else
        x1 = (x1+x2)/2; %replace as new start index
    end
    k1=k1+1;
end

rcu = fx;
P   = x2;
P1 = P1_tmp;
nSlot = nSlot_tmp;
end
end

%%
function [rcu, P1,nSlot] = golden_search_P1_SAMPR_KaKnown(f, START_INT, END_INT, TOL, E_Ka)
% Minimize min_L f(P,P1,L) over P1 for given P

P = END_INT;
iter= 20;                       % maximum number of iterations
tau=double((sqrt(5)-1)/2);      % golden proportion coefficient, around 0.618
k2=0;                            % number of iterations
x1=START_INT+(1-tau)*(END_INT-START_INT);             % computing x values
x2=START_INT+tau*(END_INT-START_INT);

TOL_nSlot = 1;
nSlot_l = E_Ka;
nSlot_u = 6*E_Ka;

[f_x1,nSlot_1] = golden_search_nSlot_SAMPR_KaKnown(f,P,x1,nSlot_l,nSlot_u,TOL_nSlot); %f(P,x1); %     % computing values in x points
[f_x2,nSlot_2] = golden_search_nSlot_SAMPR_KaKnown(f,P,x2,nSlot_l,nSlot_u,TOL_nSlot); %f(P,x2);

while ((abs(END_INT-START_INT)>TOL) && (k2<iter))
    if f_x1 < f_x2 %function smaller in x1 than in x2
        END_INT = x2; %make interval smaller
        x2 = x1; %set new end of interval
        x1 = START_INT+(1-tau)*(END_INT-START_INT); %find new beginning
        f_x2 = f_x1;%already have value in x1
        nSlot_2 = nSlot_1;
%         [fMD_x1,fFA_x1] = f(P,x1,2*E_Ka);
        [f_x1,nSlot_1] = ...
            golden_search_nSlot_SAMPR_KaKnown(f, P,x1,nSlot_l,nSlot_u,TOL_nSlot); %f(P,x1); %%compute new value for new beginning
    else
        START_INT=x1; %function smaller in x2 than in x1 so set new start indx to x1
        x1 = x2; %replace as new start index
        x2=START_INT+tau*(END_INT-START_INT); %compute new end index
        f_x1= f_x2;
        nSlot_1 = nSlot_2;
%         [fMD_x2,fFA_x2] = f(P,x2,2*E_Ka);
        [f_x2,nSlot_2] = ...
            golden_search_nSlot_SAMPR_KaKnown(f, P,x2,nSlot_l,nSlot_u,TOL_nSlot); %f(P,x2);
    end
    k2=k2+1;
end

if f_x1 < f_x2
    P1=x1;
    rcu = f_x1;
    nSlot = nSlot_1;
else
    P1=x2;
    rcu = f_x2;
    nSlot = nSlot_2;
end

end

%%
function [rcu, nSlot] = golden_search_nSlot_SAMPR_KaKnown(f, P,P1,START_INT,END_INT,TOL)
% Minimize f(P,P1,L) over L for given P and P1

iter= 20;                       % maximum number of iterations
tau=double((sqrt(5)-1)/2);      % golden proportion coefficient, around 0.618
k2=0;                            % number of iterations
x1=max(round(START_INT+(1-tau)*(END_INT-START_INT)),0);             % computing x values
x2=round(START_INT+tau*(END_INT-START_INT));

f_x1 = f(P,P1,x1);    % computing values in x points
f_x2 = f(P,P1,x2); 

while ((abs(END_INT-START_INT)>TOL) && (k2<iter))
    if f_x1 < f_x2 %function smaller in x1 than in x2
        END_INT = x2; %make interval smaller
        x2 = x1; %set new end of interval
        x1 = max(round(START_INT+(1-tau)*(END_INT-START_INT)),0); %find new beginning
        f_x2 = f_x1;%already have value in x1
        f_x1 = f(P,P1,x1); %compute new value for new beginning
    else
        START_INT=x1; %function smaller in x2 than in x1 so set new start indx to x1
        x1 = x2; %replace as new start index
        x2 = round(START_INT+tau*(END_INT-START_INT)); %compute new end index
        f_x1= f_x2;
        f_x2 = f(P,P1,x2); % 
    end
    k2=k2+1;
end

if f_x1 < f_x2
    nSlot=x1;
    rcu = f_x1;
else
    nSlot=x2;
    rcu = f_x2;
end
end