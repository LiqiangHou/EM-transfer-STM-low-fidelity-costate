%--------------------------------------------------------------------------
function [r,v] = conversion (E)
%
% Parametri orbitali

muSUN=  1.327124280000000e+011;       %gravitational parameter for the sun

a=E(1);
e=E(2);
i=E(3);
omg=E(4);
omp=E(5);
EA=E(6);


% Grandezze definite nel piano dell'orbita

b=a*sqrt(1-e^2);
n=sqrt(muSUN/a^3);

xper=a*(cos(EA)-e);
yper=b*sin(EA);

xdotper=-(a*n*sin(EA))/(1-e*cos(EA));
ydotper=(b*n*cos(EA))/(1-e*cos(EA));

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

r=R*[xper;yper;0];
v=R*[xdotper;ydotper;0];







