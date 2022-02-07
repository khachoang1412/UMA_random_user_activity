function [eps_MD, eps_FA, P,P1] = binary_search_P_MDFA(f,x1,x2,TARGET_MD,TARGET_FA,TOL,fixP1)
% Search the value of P such that eps_MD \in [TARGET_MD - TOL, TARGET_MD]
% and eps_FA \in [TARGET_FA - TOL, TARGET_FA] where eps_MD and eps_FA are 
% obtained from a f(P,P1), which is the RCU bound in Theorem 1 of 
%
% [1] K.-H. Ngo, A. Lancho, G. Durisi and A. G. i. Amat, "Unsourced
% Multiple Access With Random User Activity," submitted to IEEE Trans. Inf.
% Theory, Jan. 2022.

iter= 20;                       % maximum number of iterations
k1=0;                            % number of iterations

%% We can either optimize over P1 or choose a fixed value of P1. The latter 
%% is faster. Our experiments suggest that the optimal P1 is around 0.96*P
fracP1 = 0.96;
if fixP1 == 0
    if TARGET_MD == TARGET_FA
        search_P1 = @(x) golden_search_P1_MDFA(f,1e-9,x,x/100,'max',[],[]);
    else
        % Calculate the weights based on the target MD and FA probabilities 
        weight_MD = 1/(1 + TARGET_MD/TARGET_FA);
        weight_FA = 1/(1 + TARGET_FA/TARGET_MD);
        
        search_P1 = @(x) golden_search_P1_MDFA(f,1e-9,x,x/100,'weighted',weight_MD, weight_FA);
    end
else
    search_P1 = @(x) compute_MDFA_fixedP1(f, x, fracP1);
end

%% Search for P
[eps_MD_tmp,eps_FA_tmp,P1_tmp] = search_P1(x2);
[eps_MD_tmp1,eps_FA_tmp1,P1_tmp1] = search_P1(x1);
if x1 == x2
    P = x1;
    eps_MD = eps_MD_tmp;
    eps_FA = eps_FA_tmp;
    P1 = P1_tmp;
elseif eps_MD_tmp > TARGET_MD || eps_FA_tmp > TARGET_FA
    warning('Impossible to achieve the target within the given range of power :( ');
    P = x2;
    eps_MD = eps_MD_tmp;
    eps_FA = eps_FA_tmp;
    P1 = P1_tmp;
elseif eps_MD_tmp1 < TARGET_MD || eps_FA_tmp1 < TARGET_FA
    warning('All power values in the range can achieve the target :) ');
    P = x1;
    eps_MD = eps_MD_tmp1;
    eps_FA = eps_FA_tmp1;
    P1 = P1_tmp1;
else
    [fx_MD,fx_FA,P1_tmp] = search_P1((x1+x2)/2); 

    while ~((TARGET_MD >= fx_MD && fx_MD >= TARGET_MD - TOL && TARGET_FA >= fx_FA) || ...
            (TARGET_FA >= fx_FA && fx_FA >= TARGET_FA - TOL && TARGET_MD >= fx_MD) ...
            || (k1>iter)) 
        if k1 > 0
            [fx_MD,fx_FA,P1_tmp] = search_P1((x1+x2)/2); 
        end
        if TARGET_MD > fx_MD && TARGET_FA > fx_FA
            x2 = (x1+x2)/2; % set new end of interval        
        else
            x1 = (x1+x2)/2; % replace as new start index
        end
        k1 = k1+1;
    end

    eps_MD = fx_MD;
    eps_FA = fx_FA;
    P   = x2;
    P1 = P1_tmp;
end
end

function [eps_MD, eps_FA, P1] = compute_MDFA_fixedP1(f, P, fracP1)
    P1 = P*fracP1;
    [eps_MD, eps_FA] = f(P,P1);
end