function main_STM_Earth_Mars_low_fidelity_transfer
clear all;

format short g;
% 

%
% ------------
% lb = [-1.0*ones(1,6) 0.0];
% ub = [ 1.0*ones(1,6) 1.0];
% lambda_0 = 0.1 *ones(1,7);
% nonlcon = [];
% ------------
% low fidelity

% feq_data =  -5.7719e-07   2.8104e-06  -3.0679e-07   8.5199e-07  -3.1958e-06  -1.4299e-07   5.3608e-08       474.43       1341.9   4.3918e-06
% x =     -0.18122     -0.74713      0.44599     -0.60143      0.36878      0.17906    0.0019165     0.038467     0.049746    0.0027327
% fval =   4.3918e-06
%    73        1249    4.391786e-06     0.000e+00     2.929e-11     1.144e-15     5.472e-01  

% ode113 may lead some difference. unrecommended
   %---------------------------------------------
% high-fidelity solution

% feq_data =    7.4847e-06  -4.1743e-06  -7.2155e-06  -1.5277e-05  -6.0892e-06  -4.6833e-06   1.7962e-05       474.43       1350.9   2.7213e-05
% 
% x =   -0.28526     -0.66883      0.18613     -0.55339      0.45125      0.17403     0.010007     0.029477     0.030577   -0.0012208
% fval =
%    2.7213e-05
%    75        1306    2.721268e-05     0.000e+00     1.636e-14     1.431e-15     2.549e+00  



%---------------------------------------------------
% low fiedlity  
% % % eq.2.61, the switching function
% [u,alpha]     = Hamiltonian_switch(x,lambda_0,const);
% u = 1;
% lb = [-1.0*ones(1,7) ];
% ub = [ 1.0*ones(1,7) ];
% lambda_0 = 0.2 *ones(1,7);
%
% feq_data =
%   -4.2332e-07   2.1726e-06    4.824e-07   6.7471e-07  -3.8057e-06  -1.0518e-06   4.8284e-07       474.43       1342.5   4.6271e-06
% 
%   107        1617    4.626973e-06     0.000e+00     4.184e-11     1.513e-15     1.841e+00  
% x =
%      -0.14019     -0.75479     -0.01694     -0.48182      0.47409      0.17477    0.0058159     0.037689     0.048614    0.0055288
% 
% fval =
% 
%     4.627e-06




% initial trial of lambda 


global epsl
global epsl_t

% the optimal

addpath('asteroid propagation');
addpath('convolution-branch-nov-27');


TOF_days =  474.43;

% const values
const          = initial_const_transport(TOF_days);


%-----
epsl = 5.0e-4;
epsl_t = 1.0e-5;
%-----




%
fmincon_conv_pde_low_thrust(const);

return
end


% 
function fmincon_conv_pde_low_thrust(const)



% % Available algorithms: 'interior-point', 'sqp', 'sqp-legacy', 'active-set', and 'trust-region-reflective'.
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','TolFun',1.0e-12,'TolX',1.0e-12,'StepTolerance',1.0e-15,'PlotFcns',@optimplotfval,'MaxFunctionEvaluations',2000);
% options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','TolFun',1.0e-12,'TolX',1.0e-12,'StepTolerance',1.0e-15,'PlotFcns',@optimplotfval,'MaxFunctionEvaluations',2000);

% 
A = [];
b = [];
Aeq = [];
beq = [];
%---------------
lb = [-1.0*ones(1,7) ];
ub = [ 1.0*ones(1,7) ];

% ------------

% initial trial of lambda 

lambda_0 = 0.2 *ones(1,7);
nonlcon = [];

x0      = [lambda_0];

% paramteres of the initial velocity


vin_alpha_delta_0  =  [ 0.0,  0.0,         0.0  ];
vin_alpha_delta_lb =  [-0.05,  deg2rad(-18), deg2rad(-5)];
vin_alpha_delta_ub =  [ 0.05,  deg2rad( 18), deg2rad( 5)];



% -----

%


x0 = [x0,vin_alpha_delta_0];
lb = [lb,vin_alpha_delta_lb];
ub = [ub,vin_alpha_delta_ub];


%
fun =  @(lambda)obj_fun(lambda,const)

[x,fval]       = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)




return
end

function [c,ceq] = mycon(x)

%  costate
lambda = x(1:7);
norm_lambda = norm(lambda);

% constraints
c   = [];...                 % Compute nonlinear inequalities at x.
ceq = norm_lambda - 1.0; ...   % Compute nonlinear equalities at x.

return
end


function x0 = Initial_state_costate_state(const,vin_alpha_delta)



%-----

delta_v0  = 0.20;          % AU unit, not regular km/s
alpha0 = deg2rad(140);
%
delta0 = deg2rad(3) ;

%-----

v0_inf = delta_v0 + vin_alpha_delta(1);                                
alpha0 = alpha0   + vin_alpha_delta(2);
delta0 = delta0   + vin_alpha_delta(3);

% boundary condition, sec.2. 


[r0_earth,v0_earth,E0_earth] = pleph_an (const.mjd2000, 3);

% sate initial condition
u_alpha_delta = [...
    cos(alpha0)*cos(delta0) 
    sin(alpha0)*cos(delta0) 
    sin(delta0)              ];


x0 = [r0_earth/const.l_ref
      v0_earth/const.v_ref + v0_inf*u_alpha_delta
      const.m0/const.m_ref  ];                         % m0




% For computing convenience, the variables are adimensionalized by
% astronomycal unit (AU, 149597870.66 km ) and spacecraft initial mass
% $m_0$. Through the definition of the reference velocity as
% $v_ref}=\sqrt{\mu / AU}$, the reference time is given by the ratio
% $t_{ref}=v_{ref} / l_{ref}$.


return
end


%


function  J = obj_fun(lambda,const)

% final condition using initial value of costate
[Mu]  = transport_pde_low_thrust(lambda,const);


% obj
J  = norm(Mu);


return
end


function feq = transport_pde_low_thrust(xzero,const)
% 
vin_alpha_delta = xzero(8:10);

% initial state and costate
x0        = Initial_state_costate_state(const,vin_alpha_delta);
const.x0  = x0;

% initial state and costate

lambda_0  = xzero(1:7)'; 

%-------
y0 = [...
    x0 
    lambda_0 ];
% initila dot_x0 and u_m0

dot_x0       = initial_ode_transport_dyn(0,y0,const);
const.dot_x0 = dot_x0;




% integrator paramster
atol = 1.0e-12;
rtol = 1.0e-11;

%-------



% minstep  = 1.0e-20;
minstep  = [];
events   = [];


options = odeset('RelTol',rtol,'AbsTol',atol,'Events',events,'MinStep', minstep);
fun   = @(t,x)ode_transport_dyn(t,x,const);


% ode 
TOF_secs = const.TOF_days*86400;
tspan = [0 TOF_secs/const.t_ref];

[t,y] = ode113(fun, tspan, y0,  options);
% [t,y] = ode45(fun, tspan, y0,  options);


% boundary conditions, rendezvous at mars. sec.2.2 

tf = t(end);
yf = y(end,:);

rf = yf(1:3);
vf = yf(4:6);
mf = yf(7);


tf_days   = tf*const.t_ref/86400;

% final costate
lambda_f   = yf(8:14);
lambda_mf = lambda_f(7);


% final state

rf_mars = const.rf_mars;
vf_mars = const.vf_mars;

delta_rf = (rf - rf_mars');
delta_vf = (vf - vf_mars');

% obj vector
feq = [delta_rf, delta_vf,lambda_mf]; 

feq_data = [delta_rf, delta_vf,lambda_mf,tf_days,mf*const.m_ref,norm(feq)] 


return
end


%----
function [value,isterminal,direction] = RfEvents(t,y,const)


rf_mars = const.rf_mars;

value = y(2) - rf_mars(2);     % Detect height = 0
isterminal = 1;   % Stop the integration
direction =  1;     % The zero can be approached from either direction


return
end
%----


function dot_x = initial_ode_transport_dyn(t,y0,const)


% const value 
mu   = const.mu;

% thrust and m_dot
c_1  = const.c_1;
c_2  = const.c_2;

% state
x = y0(1:7);

r = x(1:3);
v = x(4:6);
m = x(7);

% costate, normlize to [0,1]

lambda_0 = y0(8:14);


%eq. 2.5
% -------------------------
norm_r = norm(r);
g_r    = -mu/norm_r^3*r;
% -------------------------


% % eq.2.61, the switching function
[u,alpha]     = Hamiltonian_switch(x,lambda_0,const);

u = 1;

% dot_x, eq. (3.48)

dot_rdot      =  v;
dot_vdot      =  g_r - c_1 *u/m* alpha; 
dot_m         = -1.0*c_2 *u; 


dot_x = [
     dot_rdot
     dot_vdot
     dot_m
];



return
end



%------
function const = initial_const_transport(TOF_days)


% 
% The initial mass of  the spacecraft is 1500 kg , while the propulsion system is able to
% provide a maximum thrust of 0.33 N with a $I_sp=3800 s$.

T    = 0.33;
g0   = 9.80665;
Isp  = 3800;             % sec
AU   = 149597870.66;     % km 

const.mjd2000   = 4260.62 ;                           % MJD2000
const.m0        = 1500;


% For computing convenience, the variables are adimensionalized by
% astronomycal unit (AU, 149597870.66 km ) and spacecraft initial mass
% $m_0$. Through the definition of the reference velocity as
% $v_ref}=\sqrt{\mu / AU}$, the reference time is given by the ratio
% $t_{ref}=v_{ref} / l_{ref}$.
muSUN=1.32712428e+11;       %Gravitational constant of Sun


const.l_ref  = AU;
const.v_ref  = sqrt(muSUN / AU);
const.t_ref  = const.l_ref / const.v_ref;
const.m_ref  = const.m0;
const.T_ref  = 1000*const.m_ref*(muSUN / AU^2);        % kg*km/s^2 -> kg*m/s^2


%
% eq.2.5
c_1  = T;
c_2  = T/(Isp * g0);

% adimensionalize 
const.c_1  = c_1 / const.T_ref;
const.c_2  = c_2 /(const.m_ref / const.t_ref);
const.mu   = 1.0;


% time of flight
const.TOF_days = TOF_days;


% target ephemeris 

[rf_mars, vf_mars, Ef_mars]  = pleph_an ((const.mjd2000 + TOF_days), 4);

% canonical unit
const.rf_mars = rf_mars/const.l_ref ;
const.vf_mars = vf_mars/const.v_ref ;
mf_target = 1339.5;
const.xf = [const.rf_mars;const.vf_mars; mf_target/const.m_ref];

%
const. lambda_mf =0.0;
return
end


function dot_y = ode_transport_dyn(t,y,const)

% const value 
mu   = const.mu;

% thrust and m_dot
c_1  = const.c_1;
c_2  = const.c_2;

% state
x = y(1:7);

r = x(1:3);
v = x(4:6);
m = x(7);

% costate, should not be normlized into [0,1]!!!

lambda = y(8:14);

%---
% %eq. 2.5
norm_r = norm(r);
g_r    = -mu/norm_r^3*r;
%------

% % eq.2.61, the switching function
[u,alpha]     = Hamiltonian_switch(x,lambda,const);



% dot_x, eq. (3.48)

dot_rdot      =  v;
dot_vdot      =  g_r - c_1 *u/m* alpha; 
dot_m         = -1.0*c_2 *u; 


dot_x = [
     dot_rdot
     dot_vdot
     dot_m
];



% % ------

% dot_lambda, using the state transition

dot_lambda = dot_costate_sensitive(t,x,dot_x,lambda,const);

% output

dot_y = [
    dot_x
    dot_lambda
    ];
return
end




%



function [s,alpha] = Hamiltonian_switch(x,lambda,const)
% % eq.2.5
% c_1  = T;
% c_2  = T/(Isp * g0);

m = x(7);

% thrust and m_dot
c_1  = const.c_1;
c_2  = const.c_2;


% state and costate
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

% direction of thrust
alpha = lambda_v/norm(lambda_v);
return
end





function dot_lambda = dot_costate_sensitive(t,x,dot_x,lambda,const)
global epsl
global epsl_t

x0     = const.x0;
dot_x0 = const.dot_x0;




% The gain matrix
A = eye(7,7)*epsl; 


% 
% %----

% % ---
if(t < epsl_t)
    K      = 0.0 ;
else

    delta_x     = (x - x0);
    delta_dot_x = (dot_x - dot_x0);


    K  = -1*inv(delta_x*delta_x' + A) * (delta_x*delta_dot_x'  + 0);
end


%-----

dot_lambda = K*lambda ;



return
end

%


% 



