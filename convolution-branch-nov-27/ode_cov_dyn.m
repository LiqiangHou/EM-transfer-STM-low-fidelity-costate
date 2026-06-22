
function dot_x = ode_cov_dyn(t,x,const)
% all the vectors in the computations are set to column vectors

global  x_iminus SIGMA_iminus lambda_iminus 
global  SIGMA_lambda 
% one-step predict of the state
x_i_iminus      = x;


dot_x_i_iminus = ode_motion(x_iminus,lambda_iminus,const);
% SIGAM_i_iminus = dot_x_i_iminus'*SIGMA_iminus*dot_x_i_iminus;


% obseration of H
H_i_iminus = lambda_iminus'*dot_x_i_iminus;

% one-step precit of H
dot_x_i = ode_motion(x_i_iminus,lambda_iminus,const);
H_i     = lambda_iminus'*dot_x_i;


% covaraiance of H
inv_SIGMA = pinv(SIGMA_iminus, 1.0e-10);
H_x       = partial_H_x(dot_x_i_iminus,inv_SIGMA,x_iminus,x_i_iminus);

SIGMA_H  = H_x * SIGMA_iminus * H_x';
SIGMA_XH = H_x * SIGMA_iminus;
SIGMA_HX = SIGMA_XH';


% covraiance of lambda
lambda_x      = partial_lambda_x(dot_x_i_imnius,inv_SIGMA,x_iminus,x_i_iminus);
SIGMA_lambda  = lambda_x * SIGMA_iminus;


%
inv_SIGMA_H = pinv((SIGMA_H + epsl), 1.0e-10);

% conditiona predit and covaraince 
x_i     = x_i_iminus     + SIGMA_XH*inv_SIGMA_H*(H_i - H_i_iminus); 
SIGAM_i = SIGAM_i_iminus - SIGMA_XH*inv_SIGMA_H*SIGMA_HX; 


% output
lambda_i = costate_and_likelihood(dot_x_i_iminus,x_i_imnius,x_i,SIGMA);
dot_x    = ode_motion(x_i,lambda_i,const);


% update
lambda_iminus = lambda_i;
x_iminus      = x_i;
SIGMA_iminus  = SIGAM_i;

return
end


%%


function dot_x = ode_motion(x,lambda,const)
% const value 
mu   = const.mu;


% thrust and m_dot
c_1  = const.c_1;
c_2  = const.c_2;

r = x(1:3);
v = x(4:6);
m = x(7);
%

%eq. 2.5
norm_r = norm(r);
g_r    = -mu/norm_r^3*r;

% switching function and thuts direction
[s,alpha] = Hamiltonian_switch(x,lambda,const);

% eq. (3.48)

dot_rdot      =  v;
dot_vdot      =  g_r - s*c_1 /m*alpha; 
dot_m         = -1.0*c_2 *s; 

%

dot_x = [
     dot_rdot
     dot_vdot
     dot_m ];
return
end





function lambda = costate_and_likelihood(dot_x,x_i_imnius,x_i_bar,SIGMA)

%  
inv_SIGMA = pinv(SIGMA, 1.0e-10);

% likleihood
L_X = 0.5 * (x_i_imnius - x_i_bar)*inv_SIGMA*(x_i_imnius - x_i_bar)';

% costate
lambda = dot_x*L_X;

return
end




%%


function lambda_x = partial_lambda_x(dot_x,inv_SIGMA,x_i_iminus,x_i)

% partial of lmbda w.r.t x, n\timesn matrix
lambda_x = dot_x' * (inv_SIGMA * (x_i_iminus - x_i));

return
end


function H_x = partial_H_x(dot_x,inv_SIGMA,x_i_iminus,x_i)
lambda_x = partial_lambda_x(dot_x,inv_SIGMA,x_i_iminus,x_i);

% partial of H w.r.t x
H_x = lambda_x' * dot_x;

% row vector
H_x = H_x';

return
end
