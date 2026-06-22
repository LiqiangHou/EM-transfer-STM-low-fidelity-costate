

function [s,alpha] = Hamiltonian_switch(x,lambda,const)
% % eq.2.5
% c_1  = T;
% c_2  = T/(Isp * g0);

m = x(7);


% thrust and m_dot
c_1  = const.c_1;
c_2  = const.c_2;


% state and costate
lambda_r = lambda(1 : 3);
lambda_v = lambda(4 : 6);
lambda_m = lambda(7);


% u\in[0,1]

% eq.2.61

% Defining the switching function as

rho = 1 - c_1/(c_2 * m)*norm(lambda_v) - lambda_m;

if(rho > 0)
    s = 0.0;
else
    s = 1.0;
end

% thrust direction
alpha = lambda_v/norm(lambda_v);

return
end