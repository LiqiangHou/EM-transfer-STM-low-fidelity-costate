
function J = obj_low_thrust(X_ZERO)
global const


lambda_0 = X_ZERO(1:7);
TOF_days = X_ZERO(8)



% initial value of the state and costate
x0 = Initial_state_costate(const); 
x0 = [x0
      lambda_0'];


%ode system function
TOF_secs  =  TOF_days*86400;

tspan = [0 TOF_secs/const.t_ref];


opts  = odeset(RelTol=1e-12, AbsTol=1e-12,Events=@myEventsFcn);
% opts  = odeset(RelTol=1e-12, AbsTol=1e-12);
fun   = @(t,x)ode_dyn(t,x,const);

[t,x] = ode45(fun, tspan, x0,  opts);




% boundary conditions, rendezvous at mars. sec.2.2 
tf = t(end);
tf_days   = tf*const.t_ref/86400

xf = x(end,:);
rf = xf(1:3)*const.l_ref;
vf = xf(4:6)*const.v_ref;
mf = xf(7)  *const.m_ref
lambda_mf = xf(14)
m_f      = 1295.04;
mjdf     = const.mjd2000 + TOF_days;
[rf_mars, vf_mars, Ef_mars]  = pleph_an (mjdf, 4);


% TOF_days
[rf;rf_mars';vf;vf_mars'];
feq = [(rf-rf_mars')/const.l_ref, (vf-vf_mars')/const.v_ref,lambda_mf] 
delta_rf = norm(rf - rf_mars');
delta_vf = norm(vf - vf_mars');
J = delta_rf/const.l_ref + delta_vf/const.v_ref*1.6 + abs(lambda_mf); 





return
end




function [position,isterminal,direction] = myEventsFcn(t,x)
global const

rf_mars = const.rf_mars;
vf_mars = const.vf_mars;

 % The value that we want to be zero

position = [ x(2)*const.l_ref - rf_mars(2)];
isterminal = 1;  % Halt integration 
direction =  1;   % The zero can be approached from either direction
end

