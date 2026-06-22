
%--------------------------------------------------------------------------
function E=M2E(M,e)
%
i=0;
tol=1e-10;
err=1;
E=M+e*cos(M);   %initial guess
while err>tol && i<100
    i=i+1;
    Enew=E-(E-e*sin(E)-M)/(1-e*cos(E));
    err=abs(E-Enew);
    E=Enew;
end
