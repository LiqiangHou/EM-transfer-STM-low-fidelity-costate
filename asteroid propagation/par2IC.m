

%--------------------------------------------------------------------------
function [r0,v0]=par2IC(E,mu)
%
%Usage: [r0,v0] = IC2par(E,mu)
%
%Outputs:
%           r0:    column vector for the position
%           v0:    column vector for the velocity
%
%Inputs:
%           E:     Column Vectors containing the six keplerian parameters,
%                  (a (negative fr hyperbolas),e,i,OM,om,Eccentric Anomaly 
%                    or Gudermannian if e>1)
%           mu:    gravitational constant
%
%Comments:  The parameters returned are, of course, referred to the same
%ref. frame in which r0,v0 are given. a can be given either in kms or AUs,
%but has to be consistent with mu.All the angles must be given in radians. 

a=E(1);
e=E(2);
i=E(3);
omg=E(4);
omp=E(5);
EA=E(6);


% Grandezze definite nel piano dell'orbita
if e<1
    b=a*sqrt(1-e^2);
    n=sqrt(mu/a^3);

    xper=a*(cos(EA)-e);
    yper=b*sin(EA);

    xdotper=-(a*n*sin(EA))/(1-e*cos(EA));
    ydotper=(b*n*cos(EA))/(1-e*cos(EA));
else
    b=-a*sqrt(e^2-1);
    n=sqrt(-mu/a^3);   
    dNdzeta=e*(1+tan(EA)^2)-(1/2+1/2*tan(1/2*EA+1/4*pi)^2)/tan(1/2*EA+1/4*pi);
    
    xper = a/cos(EA)-a*e;
    yper = b*tan(EA);
    
    xdotper = a*tan(EA)/cos(EA)*n/dNdzeta;
    ydotper = b/cos(EA)^2*n/dNdzeta;
end

% Matrice di trasformazione da perifocale a ECI

R(1,1)=cos(omg)*cos(omp)-sin(omg)*sin(omp)*cos(i);
R(1,2)=-cos(omg)*sin(omp)-sin(omg)*cos(omp)*cos(i);
R(1,3)=sin(omg)*sin(i);
R(2,1)=sin(omg)*cos(omp)+cos(omg)*sin(omp)*cos(i);
R(2,2)=-sin(omg)*sin(omp)+cos(omg)*cos(omp)*cos(i);
R(2,3)=-cos(omg)*sin(i);
R(3,1)=sin(omp)*sin(i);
R(3,2)=cos(omp)*sin(i);
R(3,3)=cos(i);

% Posizione nel sistema inerziale 

r0=R*[xper;yper;0];
v0=R*[xdotper;ydotper;0];

