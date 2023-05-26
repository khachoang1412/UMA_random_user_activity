function [eps_MD, eps_FA, P1] = golden_search_P1_MDFA(f, START_INT, END_INT, TOL, type, weight_MD, weight_FA)
% Minimize a weighted sum of eps_MD and eps_FA over P1 \in [0,P], where 
% eps_MD and eps_FA are given by f(P,P1), which is the RCU bound in Theorem
% 1 of 
%
% [1] K.-H. Ngo, A. Lancho, G. Durisi and A. G. i. Amat, "Unsourced
% Multiple Access With Random User Activity," submitted to IEEE Trans. Inf.
% Theory, Jan. 2022.

if strcmp(type,'max')
    g = @(fMD,fFA) max(fMD,fFA);
elseif strcmp(type,'weighted')
    g = @(fMD,fFA) fMD*weight_MD + fFA*weight_FA;
end

P = END_INT;                  % power budget

tau = double((sqrt(5)-1)/2);  % golden proportion coefficient, around 0.618

iter = 20;                    % maximum number of iterations
k2 = 0;                       % iteration index

x1 = START_INT+(1-tau)*(END_INT-START_INT);            
x2 = START_INT+tau*(END_INT-START_INT);

[fMD_x1,fFA_x1] = f(P,x1); 
[fMD_x2,fFA_x2] = f(P,x2);

while ((abs(END_INT-START_INT)>TOL) && (k2<iter))
    if g(fMD_x1, fFA_x1) < g(fMD_x2, fFA_x2) %function smaller in x1 than in x2
        END_INT = x2; %make interval smaller
        x2 = x1; %set new end of interval
        x1 = START_INT+(1-tau)*(END_INT-START_INT); %find new beginning
        fMD_x2 = fMD_x1;%already have value in x1
        fFA_x2 = fFA_x1;
        [fMD_x1,fFA_x1] = f(P,x1); %%compute new value for new beginning
    else
        START_INT=x1; %function smaller in x2 than in x1 so set new start indx to x1
        x1 = x2; %replace as new start index
        x2=START_INT+tau*(END_INT-START_INT); %compute new end index
        fMD_x1= fMD_x2;
        fFA_x1 = fFA_x2;        
        [fMD_x2,fFA_x2] = f(P,x2);
    end
    k2=k2+1;
end

% Return the solution
if g(fMD_x1, fFA_x1) < g(fMD_x2, fFA_x2)
    P1 = x1;
    eps_MD = fMD_x1;
    eps_FA = fFA_x1;
else
    P1 = x2;
    eps_MD = fMD_x2;
    eps_FA = fFA_x2;
end

end
