function [epsilon, P1] = golden_search(f, START_INT, END_INT, TOL)
% minimize f(P1) over P1. epsilon = min_{P1} f(P1)

P = END_INT;
iter = 20;                      % maximum number of iterations
tau = double((sqrt(5)-1)/2);    % golden proportion coefficient, around 0.618
k2 = 0;                         % iteration index

x1 = START_INT+(1-tau)*(END_INT-START_INT);             
x2 = START_INT+tau*(END_INT-START_INT);

f_x1 = f(P,x1);                    
f_x2 = f(P,x2);

while ((abs(END_INT-START_INT)>TOL) && (k2<iter))
    if f_x1 < f_x2 %function smaller in x1 than in x2
        END_INT = x2; % make interval smaller
        x2 = x1; % set new end of interval
        x1 = START_INT+(1-tau)*(END_INT-START_INT); % find new beginning
        f_x2 = f_x1; % already have value in x1
        f_x1 = f(P,x1); % compute new value for new beginning
    else
        START_INT=x1; % function smaller in x2 than in x1 so set new start indx to x1
        x1 = x2; % replace as new start index
        x2=START_INT+tau*(END_INT-START_INT); % compute new end index
        f_x1= f_x2;
        f_x2 = f(P,x2);
    end
    k2=k2+1;
end

if f_x1 < f_x2
    P1 = x1;
    epsilon = f_x1;
else
    P1 = x2;
    epsilon = f_x2;
end

end
