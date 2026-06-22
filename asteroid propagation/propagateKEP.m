%--------------------------------------------------------------------------
function [r,v] = propagateKEP(r0,v0,t,mu)
%
%Usage: [r,v] = propagateKEP(r0,v0,t)
%
%Inputs:
%           r0:    column vector for the non dimensional position
%           v0:    column vector for the non dimensional velocity
%           t:     non dimensional time
%
%Outputs:
%           r:    column vector for the non dimensional position
%           v:    column vector for the non dimensional velocity
%
%Comments:  The function works in non dimensional units, it takes an
%initial condition and it propagates it as in a kepler motion analytically.
%
%The matrix DD will be almost always the unit matrix, except for orbits
%with little inclination in which cases a rotation is performed so that
%par2IC is always defined
DD=eye(3);
h=vett(r0,v0);
ih=h/norm(h);
if abs(abs(ih(3))-1)<1e-3         %the abs is needed in cases in which the orbit is retrograde,
                                  %that would held ih=[0,0,-1]!!
    DD=[1,0,0;0,0,1;0,-1,0];      %Random rotation matrix that make the Euler angles well defined for the case
    r0=DD*r0;                     %For orbits with little inclination another ref. frame is used.
    v0=DD*v0;
end

E=IC2par(r0,v0,mu);  

M0=E2M(E(6),E(2));
if E(2)<1
    M=M0+sqrt(mu/E(1)^3)*t;
else
    M=M0+sqrt(-mu/E(1)^3)*t;
end
E(6)=M2E(M,E(2));
[r,v]=par2IC(E,mu);

r=DD'*r;                    
v=DD'*v;

