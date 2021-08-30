function a = newton(a,theta)

tol = 5e-12; %% tolerance for Newton's method

F = F_steady_states(a,theta);

nF = norm(F,1);
display(['At the beginning ||F|| = ',num2str(nF)])

k=0;

while (k<=20) && (nF > tol)
    DF = DF_steady_states(a,theta);
    a = a - DF\F;
    F = F_steady_states(a,theta);
    nF = norm(F);
    display(['||F|| = ',num2str(nF),', ||DF^(-1)|| = ',num2str(norm(inv(DF),1))])
    k = k+1;                            
end

end