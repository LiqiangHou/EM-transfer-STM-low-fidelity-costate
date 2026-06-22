
%--------------------------------------------------------------------------%
function E=IC2par(r0,v0,mu)
%
%Usage: E = IC2par(r0,v0,mu)
%
%Inputs:
%           r0:    column vector for the position
%           v0:    column vector for the velocity
%
%Outputs:
%           E:     Column Vectors containing the six keplerian parameters,
%                  (a (negative for hyperbolas),e,i,OM,om,Eccentric Anomaly
%                  (or Gudermannian whenever e>1))
%
%Comments:  The parameters returned are, of course, referred to the same
%ref. frame in which r0,v0 are given. Units have to be consistent, and
%output angles are in radians
%The algorithm used is quite common and can be found as an example in Bate,
%Mueller and White. It goes singular for zero inclination and for ni=pi
%Note also that the anomaly in output ranges from -pi to pi
%Note that a is negative for hyperbolae

k=[0,0,1]';
h=vett(r0,v0)';
p=h'*h/mu;
n=vett(k,h)';
n=n/norm(n);
R0=norm(r0);
evett=vett(v0,h)'/mu-r0/R0;
e=evett'*evett;
E(1)=p/(1-e);
E(2)=sqrt(e);
e=E(2);
E(3)=acos(h(3)/norm(h));
E(5)=(acos((n'*evett)/e));
if evett(3)<0
    E(5)=2*pi-E(5);
end
E(4)=acos(n(1));
if n(2)<0
    E(4)=2*pi-E(4);
end
ni=real(acos((evett'*r0)/e/R0)); %real is to avoid problems when ni~=pi
if (r0'*v0)<0
    ni=2*pi-ni;
end
EccAn=ni2E(ni,e);
E(6)=EccAn;


% column vector
E =E';



%--------------------------------------------------------------------------
function E=ni2E(ni,e)
if e<1
    E=2*atan(sqrt((1-e)/(1+e))*tan(ni/2)); %algebraic kepler's equation
else
    E=2*atan(sqrt((e-1)/(e+1))*tan(ni/2)); %algebraic equivalent of kepler's equation in terms of the Gudermannian
end


