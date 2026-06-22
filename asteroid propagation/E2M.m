%--------------------------------------------------------------------------
function M=E2M(E,e)
%
%Transforms the eccentric anomaly to mean anomaly. All i/o in radians  
%
%Usage: M=E2M(E,e)
%
%Inputs :   E : Eccentric anomaly or Gudermannian if e>1
%           e : Eccentricity of considered orbit
%
%Output :   M : Mean anomaly or N if e>1

if e<1 %Ellipse, E is the eccentric anomaly
    M=E-e*sin(E);
else  %Hyperbola, E is the Gudermannian
    M=e*tan(E)-log(tan(E/2+pi/4));
end
    
