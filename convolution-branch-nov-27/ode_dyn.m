

function dot_x = ode_dyn(t,x,const)
% const value 
mu   = const.mu;


% thrust and m_dot
c_1  = const.c_1;
c_2  = const.c_2;

r = x(1:3);
v = x(4:6);
m = x(7);

lambda_r = x(8  : 10);
lambda_v = x(11 : 13);
lambda_m = x(14);

% where $r}$ denotes the distance between the spacecraft and
% the primary and $\mu$ is the Sun gravitational constant $\mu= 
% 1.3271244001810^{11} km^3 /s^2)$.


%

%eq. 2.5
norm_r = norm(r);
g_r    = -mu/norm_r^3*r;


% eq.2.61, the switching function 
lambda = x(8  : 14);
u     = Hamiltonian_switch(x,lambda,const);
alpha = lambda_v/norm(lambda_v);

% eq. (3.48)

dot_rdot      =  v;
dot_vdot      =  g_r - c_1 *u/m* alpha; 
dot_m         = -1.0*c_2 *u; 
dot_lambda_r  =  norm(g_r)*lambda_v - (dot(3*mu*r , lambda_v) /norm_r^5)*r ;
dot_lambda_v  = -1.0*lambda_r ;
dot_lambda_m  = -1.0*c_1*u/m^2*norm(lambda_v);



dot_x = [
     dot_rdot
     dot_vdot
     dot_m
     dot_lambda_r
     dot_lambda_v
     dot_lambda_m ];
return
end



