%--------------------------------------------------------------------------
function t = time2distance(r0,v0,rtarget)
%
%Usage: t = time2distance(r0,v0,rtarget)
%
%Inputs:
%           r0:    column vector for the position (mu=1)
%           v0:    column vector for the velocity (mu=1)
%           rtarget: distance to be reached
%
%Outputs:
%           t:     time taken to reach a given distance
%
%Comments:  everything works in non dimensional units

r0norm = norm(r0);
if r0norm < rtarget
    out = sign(r0'*v0);
    E = IC2par(r0,v0,1);
    a = E(1); e = E(2); E0 = E(6); p = a * (1-e^2);
    %If the solution is an ellipse 
    if e<1
        ra = a * (1+e);
        if rtarget>ra
            t = NaN;
        else %we find the anomaly where the target distance is reached
            ni = acos((p/rtarget-1)/e);         %in 0-pi
            Et = ni2E(ni,e);          %in 0-pi
            if out==1
                t = a^(3/2)*(Et-e*sin(Et)-E0 +e*sin(E0));
            else
                E0 = -E0;
                t = a^(3/2)*(Et-e*sin(Et)+E0 - e*sin(E0));
            end
        end
    else %the solution is a hyperbolae
        ni = acos((p/rtarget-1)/e);         %in 0-pi
        Et = ni2E(ni,e);          %in 0-pi
        if out==1
                t = (-a)^(3/2)*(e*tan(Et)-log(tan(Et/2+pi/4))-e*tan(E0)+log(tan(E0/2+pi/4)));
            else
                E0 = -E0;
                t = (-a)^(3/2)*(e*tan(Et)-log(tan(Et/2+pi/4))+e*tan(E0)-log(tan(E0/2+pi/4)));
        end
    end
else
    t=12;
end
    