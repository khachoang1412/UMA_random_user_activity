function [epsilon, P] = binary_search(f,x1,x2,TARGET,TOL)
% search the value of P such that min_{P1} f(P,P1) is from TARGET - TOL to TARGET

iter = 20;                       % maximum number of iterations
k1 = 0;                            % iteratin index

% we use a golden search to compute min_{P1} f(P,P1)
fx = golden_search(f,1e-9,(x1+x2)/2,(x1+x2)/200); 

while (TARGET < fx || fx<TARGET - TOL) && (k1<iter || TARGET < fx)
    if k1 > 0
        fx = golden_search(f,0,(x1+x2)/2,(x1+x2)/200); 
    end
    if TARGET > fx
        x2 = (x1+x2)/2; %set new end of interval        
    else
        x1 = (x1+x2)/2; %replace as new start index
    end
    k1=k1+1;
end

epsilon = fx;
P   = x1;
end
