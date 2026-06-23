function plot_STM_Earth_Mars_low_fidelity_transfer
clear all;
close all;


format short g;





addpath('asteroid propagation');
addpath('convolution-branch-nov-27');


%-- optimal initla guess of costate and departure parameters ---
xzero =    [       -0.14019     -0.75479     -0.01694     -0.48182      0.47409      0.17477    0.0058159     0.037689     0.048614    0.0055288];

%
TOF_days =  474.43;

% const values
const          = initial_const_transport(TOF_days);

% construct the transfer
[feq,t,y,const] = transport_pde_low_thrust(xzero,const);


% 
plot_trajetcory(t,y,const);

return
end

%% ---
function plot_trajetcory(t,y,const)
t_days = t*const.t_ref/86400;

% trajectory of mars and earth
for k = 1:size(t, 1)
    mjd     = const.mjd2000 + t_days(k);


    [r_mars,  v_mars,  E_mars]   = pleph_an (mjd, 4);
    [r_earth, v_earth, E_earth]  = pleph_an (mjd, 3);


    array_r_mars(k , :)  = r_mars/const.l_ref;
    array_r_earth(k , :) = r_earth/const.l_ref;
    
end

% % font and style setting
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'DefaultAxesFontName', 'Times New Roman');
set(groot, 'DefaultAxesFontSize', 14);
set(groot, 'DefaultLineLineWidth', 1);
set(gca, 'FontSize', 14);

% start and final position
pos_earth = const.r0_earth;
pos_mars  = const.rf_mars;



% --- switch & thrust---

h_thrust = figure(1);

% switch and thrust profile

for k = 1:size(t, 1)
    [dot_y,u,alpha]      = ode_transport_dyn(t(k),y(k,:)',const);

    t_swicth_alpha(k,:)  = [t(k),u,-alpha'];
end

hold on;
plot(t_days,t_swicth_alpha(:,2), 'k-','LineWidth', 2.5);
plot(t_days,t_swicth_alpha(:,3), 'k-.','LineWidth',2.5 );
plot(t_days,t_swicth_alpha(:,4), 'k:','LineWidth', 2.5 );
plot(t_days,t_swicth_alpha(:,5), 'k--','LineWidth',2.5 );
xlabel('Time of Flight (days)');
ylabel('Thrust Profile');

axis tight;
xlim([min(t_days) max(t_days)]);
ylim padded;
legend('Switch function', 'Thrust Direction - x','Thrust Direction - y','Thrust Direction - z');

% Detect changes
[pts_switch, ~] = findchangepts(t_swicth_alpha(:,2), 'MaxNumChanges', 1);
savefig(h_thrust, 'low_fidelity_earth_mars_thrust.fig');

% --- transfer ----
h_trajector = figure(2);

hold on; 



% Earth: 
% Plot a red circle marker at the end point
plot(pos_earth(1), pos_earth(2), 'ko', ...
    'MarkerSize', 10, ...
    'MarkerFaceColor', 'k'); % 'ro' for red circle

label_text = 'Earth';
% Adjust the position for better readability, e.g., slightly above/right
text(pos_earth(1)+0.02, pos_earth(2)+0.02,label_text, ...
    'VerticalAlignment', 'bottom', ... % Align text below the point
    'HorizontalAlignment', 'left', ...  % Align text to the left of the point
    'Color', 'k', ... % 'k' for black text
    'FontSize', 14); 

% Mars: 
% Plot a red circle marker at the end point
plot(pos_mars(1), pos_mars(2), 'kp', ...
    'MarkerSize', 10, ...
    'MarkerFaceColor', 'k'); % 'bp' for red circle

label_text = 'Mars';
% Adjust the position for better readability, e.g., slightly above/right
text(pos_mars(1)+0.02, pos_mars(2)+0.02,label_text, ...
    'VerticalAlignment', 'bottom', ... % Align text below the point
    'HorizontalAlignment', 'left', ...  % Align text to the left of the point
    'Color', 'k', ... % 'k' for black text
    'FontSize', 14); 


% Plot the first and last points with a different marker (e.g., a red asterisk)
plot( 0 , 0 ,              'k+',...
array_r_mars(:,1),   array_r_mars(:,2),   'k-.',...
array_r_earth(:,1) , array_r_earth(:,2),  'k-.'...
);


% Highlight thrust segment  with a different style
plot(y(:,1) , y(:,2),       'k--', 'LineWidth', 2) ;
plot(y(pts_switch:end, 1), y(pts_switch:end , 2),'k', 'LineWidth', 4) 

legend('Erath', 'Mars','Sun','Spacecraft');


axis square; 

% % 4. Add the vector (arrow)
% quiver(y(pts_switch:end, 1), y(pts_switch:end , 2), t_swicth_alpha(pts_switch:end,3)*0.5, t_swicth_alpha(pts_switch:end,4)*0.5, 'k', 'LineWidth', 2)

savefig(h_trajector, 'low_fidelity_earth_mars_trajectory.fig');
hold off;



% Requires export_fig toolbox from File Exchange
exportgraphics(h_trajector, 'low_fidelity_earth_mars_trajectory.pdf', 'ContentType', 'vector');
exportgraphics(h_thrust, 'low_fidelity_earth_mars_thrust.pdf', 'ContentType', 'vector');


return
end


%%

function x0 = Initial_state_costate_state(const,vin_alpha_delta)

%-----
v0_inf  = 0.20 + vin_alpha_delta(1);
%-----

alpha0 = deg2rad(140);
delta0 = deg2rad(3) ;

%-----

                           
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



function [feq,t,y,const] = transport_pde_low_thrust(xzero,const)
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


% the value of u is set to 1 for initialize dot_x at departure 
u=1;

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



% initial conditon
[r0_earth,v0_earth,E0_earth] = pleph_an (const.mjd2000, 3);
const.r0_earth = r0_earth'/const.l_ref;
const.v0_earth = v0_earth'/const.l_ref;




return
end


function [dot_y,u,alpha] = ode_transport_dyn(t,y,const)

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

%---


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
%-----
epsl = 5.0e-4;
epsl_t = 1.0e-5;
%-----


x0     = const.x0;
dot_x0 = const.dot_x0;




% The gain matrix
A = eye(7,7)*epsl; 


% ---
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



