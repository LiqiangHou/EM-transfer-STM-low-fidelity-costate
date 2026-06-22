% ------------------------------------------------------------------------
% This source file is part of the 'ESA Advanced Concepts Team's			
% Space Mechanics Toolbox' software.                                       
%                                                                          
% The source files are for research use only,                              
% and are distributed WITHOUT ANY WARRANTY. Use them on your own risk.     
%                                                                                                                                                  
% Copyright (c) 2004-2007 European Space Agency
% ------------------------------------------------------------------------
% 

%Programmed by:         Dario Izzo          (ESA/ACT)
%                       Claudio Bombardelli (ESA/ACT)

%Date:                  15/03/2007
%Revision:              4
%Tested by:             D.Izzo and C.Bombardelli
%
%
%  Computes the DeltaV cost function of a Multiple Gravity Assist trajectory
%  with a Deep Space Maneuver between each planet pair
%  N.B.: All swing-bys are UNPOWERED (thrust is only present at each dsm)
%  It takes as input a sequence of planets P1,Pn and a decision vector t
%
%  The flyby/dsm sequence is:  P1/d1/P2/d2/../dn-1/Pn
%

% DECISION VECTOR DEFINITION:
%
%  t(1)=  epoch of departure MJD2000 from first planet (not necessarily earth)
%  t(2)=  magnitude hyperbolic escape velocity from first planet
%  t(3,4)= u,v variables for the hyperbolic velocity orientation wrt Earth
%  velocity at departure
%  t(5..n+3)= planet-to-planet Time of Flight [ToF] (days)
%  t(n+4..2n-2)= fraction of ToF at which the DSM occurs
%  t(2n+3..3n) = perigee fly-by radius for planets P2..Pn-1, non-dimensional wrt planetary radii
%  t(3n+1..4n-2) = rotation gamma of the bplane-component of the swingby outgoing velocity (v_rel_out)
%  [take n_r=cross(v_rel_in,v_planet_helio) if you rotate n_r by +gamma around v_rel_in
%  you obtain the projection of v_rel_out on the b-plane]
%  Vector (Vout) around the axis of the incoming swingby velocity vector (Vin)
%

%Usage:     [J,DVvec,DVarr] = mga_dsm(t,MGADSMproblem)
%Outputs:
%           J:     Cost function = depends on the problem:
%                  orbit insertion: total DV from propulsion system (V_launcher not counted)
%                  gtoc1= (s/c final mass)*v_asteroid'*vrel_ast_sc
%                  asteroid deflection-> 1/d with d=deflecion on the earth-asteroid lineofsight
%           DVvec: vector of all DV maneuvers, including the escape C3 and
%                  the arrival DV evaluaed according to the objective function
%           DVarr: Relative velocity at the arrival planet
%
%
%Inputs:
%           t:         decision vector
%               MGADSMproblem = struct array defining the problem, i.e.
%               MGADSMproblem.sequence:  planet sequence. Example [3 3 4]= [E E M]. A
%                      negative sign represents a retrograde orbit after the fly-by
%               MGADSMproblem.objective.type: type of objective function (e.g.
%                       'orbit insertion','rndv','gtoc1')
%               MGADSMproblem.objective.rp: pericentre radius of the
%                       target orbit if type = 'orbit insertion'
%               MGADSMproblem.objective.e:  eccentricity of the target
%                       orbit if type = 'orbit insertion'
%               MGADSMproblem.bounds: decision vector upper and lower
%               bounds for the optimiser
%               MGADSMproblem.yplot:     1-> plot trajectory, 0-> don't
%
%
%*********  IMPORTANT NOTE (SINGULARITY)   ************
%
% The routine is singular when the S/C relative incoming  velocity to a planet v_rel_in
% is parallel to the heliocentric velocity of that planet.
%
% One can move the singularity elsewhere by changing the definition of the
% angle gamma.
% For example one possibility is to define gamma as follows:
% [take n_r=cross(v_rel_in,r_planet_helio) and rotate n_r by +gamma around v_rel_in
% in this case is singular for v_rel_in parallel to r_planet_helio
%
%
function [J,DVvec,DVrel] = GTOC6_mga_dsm(t,problem)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sequence= problem.sequence;  % THE PLANETs SEQUENCE


switch problem.objective.type
    case 'orbit insertion'
        rp_target= problem.objective.rp; % radius of pericentre at capture
        e_target=problem.objective.e;   % Eccentricity of the target orbit at capture
    case 'total DV orbit insertion'
        rp_target= problem.objective.rp; % radius of pericentre at capture
        e_target=problem.objective.e;   % Eccentricity of the target orbit at capture
    case 'gtoc1'
        Isp=problem.objective.Isp;
        mass = problem.objective.mass;
    case 'time to AUs'
        AUdist = problem.objective.AU;
        DVtotal = problem.objective.DVtot;
        DVonboard = problem.objective.DVonboard;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%**************************************************************************
%Definition of the gravitational constants of the various planets
%(used in the powered swing-by routine) and of the sun (used in the lambert
%solver routine)
%**************************************************************************

mu(1)=22321;          %                          Mercury
mu(2)=324860;         %Gravitational constant of Venus
mu(3)=398601.19;      %                          Earth
mu(4)=42828.3;        %                          Mars
mu(5)=126.7e6;        %                          Jupiter
mu(6)=37.93951970883e6; %                        Saturn

muSUN=1.32712428e+11; %Gravitational constant of Sun


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Definition of planetari radii
%
RPL(1)=2440; % Mercury
RPL(2)=6052; % Venus
RPL(3)=6378; % Earth
RPL(4)=3397; % Mars
RPL(5)=71492;% Jupiter
RPL(6)=60330;% Saturn


%**************************************************************************
% Decision vector definition

tdep=t(1);         % departure epoch (MJD2000)
VINF=t(2);         %Hyperbolic escape velocity (km/sec)
udir=t(3);            %Hyperbolic escape velocity var1 (non dim)
vdir=t(4);            %Hyperbolic escape velocity var2 (non dim)
N=length(sequence);

%Preallocating memory increase speed
tof = zeros(N-1,1);
alpha = zeros(N-1,1);

for i=1:N-1
    tof(i)=t(i+4); %planet-to-planet Time of Flight [ToF] (days)
    alpha(i)=t(N+i+3); %fraction of ToF at which the DSM occurs
end


%If we are optimising to reach a given distance in the shortest time ('time
%to AUs) then the decision vector needs to include also r_p and bincl at
%the last planet, otherwise not

if strcmp(problem.objective.type,'time to AUs')
    rp_non_dim=zeros(N-1,1);        %initialization gains speed
    gamma=zeros(N-1,1);
    for i=1:N-1
        rp_non_dim(i)=t(i+2*N+2); % non-dim perigee fly-by radius of planets P2..Pn (i=1 refers to the second planet)
        gamma(i)=t(3*N+i);        % rotation of the bplane-component of the swingby outgoing
        % velocity  Vector (Vout) around the axis of the incoming swingby velocity vector (Vin)
    end
else
    rp_non_dim=zeros(N-2,1);        %initialization gains speed
    gamma=zeros(N-2,1);
    for i=1:N-2
        rp_non_dim(i)=t(i+2*N+2); % non-dim perigee fly-by radius of planets P2..Pn-1 (i=1 refers to the second planet)
        gamma(i)=t(3*N+i);        % rotation of the bplane-component of the swingby outgoing
        % velocity  Vector (Vout) around the axis of the incoming swingby velocity vector (Vin)
    end
end

%**************************************************************************
%Evaluation of position and velocities of the planets
%**************************************************************************

N=length(sequence);

r= zeros(3,N);
v= zeros(3,N);

muvec=zeros(N,1);
Itime = zeros(N);
dT=zeros(1,N);

seq = abs (sequence);

T=tdep;
dT(1:N-1)=tof;
figure(3)%%%%%画出行星位置
hold on 
for i=1:N
    Itime(i)=T;
    if seq(i)<10
        [r(:,i),v(:,i)]=pleph_an( T , seq(i)); %positions and velocities of solar system planets
        muvec(i)=mu(seq(i)); %gravitational constants
    else
        [r(:,i),v(:,i)]=CUSTOMeph( mjd20002jed(T) , ...
            problem.customobject(seq(i)).epoch, ...
            problem.customobject(seq(i)).keplerian , 1); %positions and velocities of custom object
        muvec(i)=problem.customobject(seq(i)).mu; %gravitational constant of custom object

    end

    T=T+dT(i);
    
    plot3(r(1,i),r(2,i),r(3,i),'go');%%%%%画出行星位置
end
plot3(0,0,0,'ro');%%%%%画出太阳位置


if strcmp(problem.objective.type,'time to AUs')
    rp=zeros(N-1,1);        %initialization gains speed
    for i=1:N-1
        rp(i)= rp_non_dim(i)*RPL(seq(i+1)); %dimensional flyby radii (i=1 corresponds to 2nd planet)
    end
else
    rp=zeros(N-2,1);        %initialization gains speed
    for i=1:N-2
        rp(i)= rp_non_dim(i)*RPL(seq(i+1)); %dimensional flyby radii (i=1 corresponds to 2nd planet)
    end
end


%**************************************************************************
%%%% FIRST BLOCK (P1 to P2)

%Spacecraft position and velocity at departure

vtemp= cross(r(:,1),v(:,1));

iP1= v(:,1)/norm(v(:,1));
zP1= vtemp/norm(vtemp);
jP1= cross(zP1,iP1);


theta=2*pi*udir;         %See Picking a Point on a Sphere
phi=acos(2*vdir-1)-pi/2; %In this way: -pi/2<phi<pi/2 so phi can be used as out-of-plane rotation

%vinf=VINF*(-sin(theta)*iP1+cos(theta)*cos(phi)*jP1+sin(phi)*cos(theta)*zP1);
vinf=VINF*(cos(theta)*cos(phi)*iP1+sin(theta)*cos(phi)*jP1+sin(phi)*zP1);

v_sc_pl_in(:,1)=v(:,1); %Spacecraft absolute incoming velocity at P1
v_sc_pl_out(:,1)=v(:,1)+vinf; %Spacecraft absolute outgoing velocity at P1

%Days from P1 to DSM1
tDSM(1)=alpha(1)*tof(1);

%Computing S/C position and absolute incoming velocity at DSM1
[rd(:,1),v_sc_dsm_in(:,1)]=propagateKEP(r(:,1),v_sc_pl_out(:,1),tDSM(1)*24*60*60,muSUN);

%Evaluating the Lambert arc from DSM1 to P2

lw=vett(rd(:,1),r(:,2));
lw=sign(lw(3));
if lw==1
    lw=0;
else
    lw=1;
end
[v_sc_dsm_out(:,1),v_sc_pl_in(:,2)]=lambertI(rd(:,1),r(:,2),tof(1)*(1-alpha(1))*24*60*60,muSUN,lw);

%First Contribution to DV (the 1st deep space maneuver)
DV=zeros(N-1,1);
DV(1)=norm(v_sc_dsm_out(:,1)-v_sc_dsm_in(:,1));


%****************************************
% INTERMEDIATE BLOCK

tDSM=zeros(N-1,1);
for i=1:N-2

    %Evaluation of the state immediately after Pi

    v_rel_in=v_sc_pl_in(:,i+1)-v(:,i+1);

    e=1+rp(i)/muvec(i+1)*v_rel_in'*v_rel_in;

    beta_rot=2*asin(1/e);              %velocity rotation

    ix=v_rel_in/norm(v_rel_in);
    % ix=r_rel_in/norm(v_rel_in);  % activating this line and disactivating the one above
    % shifts the singularity for r_rel_in parallel to v_rel_in

    iy=vett(ix,v(:,i+1)/norm(v(:,i+1)))';
    iy=iy/norm(iy);
    iz=vett(ix,iy)';
    iVout = cos(beta_rot) * ix + cos(gamma(i))*sin(beta_rot) * iy + sin(gamma(i))*sin(beta_rot) * iz;
    v_rel_out=norm(v_rel_in)*iVout;

    v_sc_pl_out(:,i+1)=v(:,i+1)+v_rel_out;


    %Days from Pi to DSMi
    tDSM(i+1)=alpha(i+1)*tof(i+1);


    %Computing S/C position and absolute incoming velocity at DSMi
    [rd(:,i+1),v_sc_dsm_in(:,i+1)]=propagateKEP(r(:,i+1),v_sc_pl_out(:,i+1),tDSM(i+1)*24*60*60,muSUN);


    %Evaluating the Lambert arc from DSMi to Pi+1

    lw=vett(rd(:,i+1),r(:,i+2));
    lw=sign(lw(3));
    if lw==1
        lw=0;
    else
        lw=1;
    end
    [v_sc_dsm_out(:,i+1),v_sc_pl_in(:,i+2)]=lambertI(rd(:,i+1),r(:,i+2),tof(i+1)*(1-alpha(i+1))*24*60*60,muSUN,lw);

    %DV contribution
    DV(i+1)=norm(v_sc_dsm_out(:,i+1)-v_sc_dsm_in(:,i+1));

end

rd
plot3(rd(1,:),rd(2,:),rd(3,:))


%************************************************************************
% FINAL BLOCK
%
%1)Evaluation of the arrival DV
%
DVrel=norm(v(:,N)-v_sc_pl_in(:,N)); %Relative velocity at target planet


switch problem.objective.type
    case 'orbit insertion'
        DVper=sqrt(DVrel^2+2*muvec(N)/rp_target);  %Hyperbola
        DVper2=sqrt(2*muvec(N)/rp_target-muvec(N)/rp_target*(1-e_target)); %Ellipse
        DVarr=abs(DVper-DVper2);
    case 'total DV orbit insertion'
        DVper=sqrt(DVrel^2+2*muvec(N)/rp_target);  %Hyperbola
        DVper2=sqrt(2*muvec(N)/rp_target-muvec(N)/rp_target*(1-e_target)); %Ellipse
        DVarr=abs(DVper-DVper2);
    case 'rndv'
        DVarr = DVrel;
    case 'total DV rndv'
        DVarr = DVrel;
    case 'gtoc1'
        DVarr = DVrel;
    case 'time to AUs'  %no DVarr is considered
        DVarr = 0;
end

DV(N)=DVarr;


%
%**************************************************************************
%Evaluation of total DV spent by the propulsion system
%**************************************************************************

switch problem.objective.type
    case 'gtoc1'
        DVtot=sum(DV(1:N-1));
    case 'deflection demo'
        DVtot=sum(DV(1:N-1));
    otherwise
        DVtot=sum(DV);
end


DVvec=[VINF ; DV];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Finally our objective function is:
switch problem.objective.type
    case 'total DV orbit insertion'
        J= DVtot+VINF;
    case 'total DV rndv'
        J= DVtot+VINF;
    case 'orbit insertion'
        J= DVtot;
    case 'rndv'
        J= DVtot;
    case 'gtoc1'
        mass_fin = mass * exp (- DVtot/ (Isp/1000 * 9.80665));
        J = 1/(mass_fin * abs((v_sc_pl_in(:,N)-v(:,N))'* v(:,N)));
    case 'deflection demo'
        mass_fin = mass * exp (- DVtot/ (Isp/1000 * 9.80665));
        %non-dimensional units (mu=1)
        AU=149597870.66;
        VEL=sqrt(muSUN/AU);
        TIME=AU/VEL;
        %calculate the DV due to the impact
        relV = (V(:,N,1)-v(:,N));
        impactDV = (relV * mass_fin/astmass)/VEL;
        %calculate phi (see defl_radial.m)
        ir = r(:,N)/norm(r(:,N));
        iv = v(:,N)/norm(v(:,N));
        ih = vett(ir,iv)';
        ih = ih/norm(ih);
        itheta = vett(ih,ir)';
        impactDV = (impactDV'*ir) * ir + (impactDV'*itheta) * itheta; %projection on the orbital plane
        phi = acos((impactDV/norm(impactDV))'*ir);
        if (impactDV'*itheta) < 0
            phi = phi +pi;
        end
        %calculate ni0
        r0 = r(:,N)/AU;
        v0 = v(:,N)/VEL;
        E0 = IC2par( r0 , v0 , 1 );
        ni0 = E2ni(E0(6),E0(2));
        %calcuate the deflection projected in covenient coordinates
        a0 = E0(1);
        e0 = E0(2);
        M0 = E2M(E0(6),E0(2));
        numberofobs = 50;
        obstime = linspace(0,obstime,numberofobs);
        M = M0 + obstime*60*60*24/TIME*sqrt(1/a0);
        theta=zeros(numberofobs,1);
        for jj = 1:numberofobs
            theta(jj) = M2ni(M(jj),E0(2));
        end
        theta = theta-ni0;
        [dumb,dr] = defl_radial(a0,e0,ni0,phi,norm(impactDV),theta);

        [dumb,dtan] = defl_tangential(a0,e0,ni0,phi,norm(impactDV),theta);

        %calculate the deflecion on the Earth-asteroid lineofsight
        defl=zeros(3,numberofobs);
        temp=zeros(numberofobs,1);
        for i=1:numberofobs
            Tobs = T + obstime(i);
            [rast,vast]=CUSTOMeph( mjd20002jed(Tobs) , ...
                problem.customobject(seq(end)).epoch, ...
                problem.customobject(seq(end)).keplerian , 1);
            [rearth,vearth]=pleph_an( Tobs , 3);
            lineofsight=(rearth-rast)/norm((rearth-rast));
            defl(:,i) = rast / norm(rast) * dr(i) + ...
                vast / norm(vast) * dtan(i);
            temp(i) = norm(lineofsight'*(defl(:,i)));

        end

        J = 1./abs(temp)/AU;
        [J,index]=min(J);
    case 'time to AUs'
        %non dimensional units
        AU  = 149597870.66;
        V = sqrt(muSUN/AU);
        T = AU/V;
        %evaluate the state of the spacecraft after the last fly-by
        v_rel_in=v_sc_pl_in(:,N)-v(:,N);
        e=1+rp(N-1)/muvec(N)*v_rel_in'*v_rel_in;
        beta_rot=2*asin(1/e);              %velocity rotation
        ix=v_rel_in/norm(v_rel_in);
        % ix=r_rel_in/norm(v_rel_in);  % activating this line and disactivating the one above
        % shifts the singularity for r_rel_in parallel to v_rel_in
        iy=vett(ix,v(:,N)/norm(v(:,N)))';
        iy=iy/norm(iy);
        iz=vett(ix,iy)';
        iVout = cos(beta_rot) * ix + cos(gamma(N-1))*sin(beta_rot) * iy + sin(gamma(N-1))*sin(beta_rot) * iz;
        v_rel_out=norm(v_rel_in)*iVout;
        v_sc_pl_out(:,N)=v(:,N)+v_rel_out;
        t = time2distance(r(:,N)/AU,v_sc_pl_out(:,N)/V,AUdist);
        DVpen=0;
        if sum(DVvec)>DVtotal
            DVpen=DVpen+(sum(DVvec)-DVtotal);
        end
        if sum(DVvec(2:end))>DVonboard
            DVpen=DVpen+(sum(DVvec(2:end))-DVonboard);
        end

        J= (t*T/60/60/24 + sum(tof))/365.25 + DVpen*100;
        if isnan(J)
            J=100000;
        end
end








%--------------------------------------------------------------------------
%Subfunction that evaluates the time of flight as a function of x
function t=x2tof(x,s,c,lw,N)  

am=s/2;
a=am/(1-x^2);
if x<1 %ELLISSE
    beta=2*asin(sqrt((s-c)/2/a));
    if lw
        beta=-beta;
    end
    alfa=2*acos(x);
else   %IPERBOLE
    alfa=2*acosh(x);
    beta=2*asinh(sqrt((s-c)/(-2*a)));
    if lw
        beta=-beta;
    end
end
t=tofabn(a,alfa,beta,N);


%--------------------------------------------------------------------------
function t=tofabn(sigma,alfa,beta,N)
%
%subfunction that evaluates the time of flight via Lagrange expression
%
if sigma>0
    t=sigma*sqrt(sigma)*((alfa-sin(alfa))-(beta-sin(beta))+N*2*pi);
else
    t=-sigma*sqrt(-sigma)*((sinh(alfa)-alfa)-(sinh(beta)-beta));
end

%--------------------------------------------------------------------------
function v=vers(V) 
%subfunction that evaluates unit vectors
v=V/sqrt(V'*V);






