function const = initial_const()



%
% consider the case in which the spacecraft starts from
% Earth's position, with a $v_{\infty}^0=0.2 km/s$
% with respect to Earth's velocity, and arrives at Mars, with the same
% heliocentric velocity of this planet. As already explained in section 2.2
% , this value of $v_{\infty}^0$ is considered fixed, defined by the user
% or alternatively obtained from the global solution. 
% 
% The initial mass of  the spacecraft is 1500 kg , while the propulsion system is able to
% provide a maximum thrust of 0.33 N with a $I_sp=3800 s$.

T    = 0.33;
g0   = 9.80665;
Isp  = 3800;             % sec
AU   = 149597870.66;     % km 


% The search space for departure epoch and transfer times is [4000, 4300]
% MJD2000, and $[200,500]$ days respectively. The other unknowns,
% $\xi_r,\xi_v,\xi_lambda_r or\xi_lambda_v$ are sought
% in the interval $[-\infty, \infty]$. As already explained in section 3.5,
% the search space for the mass multipliers $\xi_{\lambda_m}$ is set to
% $[0,1]$. These assumptions hold for all the following test cases, except
% where otherwise specified. Finally, for this particular case, the final
% mass is sought in the interval $[500,1500]kg$.

const.mjd2000   = 4260.62 ;                           % MJD2000
const.m0        = 1500;


% For computing convenience, the variables are adimensionalized by
% astronomycal unit (AU, 149597870.66 km ) and spacecraft initial mass
% $m_0$. Through the definition of the reference velocity as
% $v_ref}=\sqrt{\mu / AU}$, the reference time is given by the ratio
% $t_{ref}=v_{ref} / l_{ref}$.
muSUN=1.32712428e+11;       %Gravitational constant of Sun


const.l_ref  = AU;
% const.t_ref  = sqrt(AU^3 / muSUN);
% const.v_ref  = const.l_ref / const.t_ref;

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
TOF_days      =  474.43;
const.TOF_days = TOF_days;

% final condition
mjdf     = const.mjd2000 + TOF_days;
[rf_mars, vf_mars, Ef_mars]  = pleph_an (mjdf, 4);

const.rf_mars = rf_mars';
const.vf_mars = vf_mars';
const.Ef_mars = Ef_mars';


const.lambda_mf = 0.0;


% initial conditon
[r0_earth,v0_earth,E0_earth] = pleph_an (const.mjd2000, 3);
const.r0_earth = r0_earth';
const.v0_earth = v0_earth';
const.E0_earth = E0_earth';


% initial conditon
[r0_earth,v0_earth,E0_earth] = pleph_an (const.mjd2000, 3);
const.r0_earth = r0_earth';
const.v0_earth = v0_earth';
const.E0_earth = E0_earth';

