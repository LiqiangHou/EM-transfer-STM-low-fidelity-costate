%%

function x0 = Initial_state_costate(const)

alpha0    = 1.230 ;
delta0    = 0.128 ;
nu0       = -3.432 ;

v0_inf    = 0.2;                                 % km/s


% boundary condition, sec.2. 


[r0_earth,v0_earth,E0_earth] = pleph_an (const.mjd2000, 3);

% sate initial condition
u_alpha_delta = [...
    cos(alpha0)*cos(delta0) 
    sin(alpha0)*cos(delta0) 
    sin(delta0)              ];


r0 = r0_earth;
v0 = v0_earth + v0_inf*u_alpha_delta;
m0 = const.m0;

x0 = [r0/const.l_ref
      v0/const.v_ref
      m0/const.m_ref  ];                         % m0
          
x0(4) = x0(4) + 0.16;

% For computing convenience, the variables are adimensionalized by
% astronomycal unit (AU, 149597870.66 km ) and spacecraft initial mass
% $m_0$. Through the definition of the reference velocity as
% $v_ref}=\sqrt{\mu / AU}$, the reference time is given by the ratio
% $t_{ref}=v_{ref} / l_{ref}$.




return
end

