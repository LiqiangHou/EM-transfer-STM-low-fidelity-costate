

function H = Hamiltonian_function(dot_x,lambda)

H = lambda*dot_x';

return
end



function H = bkup_Hamiltonian_function(dot_x,lambda,const)

% const value
mu   = const.mu;

% thrust and m_dot
c_1  = const.c_1;
c_2  = const.c_2;


% state and costate
r = x(1:3);
v = x(4:6);
m = x(7);

lambda_r = lambda(1 : 3);
lambda_v = lambda(4 : 6);
lambda_m = lambda(7);


% u\in[0,1]

% eq.2.23
norm_r =  norm(r);
alpha  = -lambda_v/norm(lambda_v);

% H = dot(lambda_r ,  v)                                    ...
%   + dot(lambda_v ,  (-mu/norm_r^3*r + s*c_1/m *alpha))    ...
%   - lambda_m* c_2 * s                                     ...
%   + c_2 * s;
return
end

