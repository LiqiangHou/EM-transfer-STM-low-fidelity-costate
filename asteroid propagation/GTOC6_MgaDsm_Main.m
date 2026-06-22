%GTOC6_mga_dsm的使用（自己编的例子）



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
%Usage:     [J,DVvec,DVarr] = GTOC6_mga_dsm(t,MGADSMproblem)
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
%                      negative sign represents a retrograde（后向） orbit after the fly-by
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

 t(1)=  60000;%epoch of departure MJD2000 from first planet (not necessarily earth);;% departure epoch (MJD2000) 
 t(2)=  3.4;%magnitude hyperbolic escape velocity from first planet;;%Hyperbolic escape velocity (km/sec)
 %t(3,4)=[1 1];% u,v variables for the hyperbolic velocity orientation wrt Earth velocity at departure
 t(3)=0;%Hyperbolic escape velocity var1 (non dim)
 t(4)=1;%Hyperbolic escape velocity var2 (non dim)
 %n=3
 %t(5,6)=[20 30];% planet-to-planet Time of Flight [ToF] (days)
 t(5)=20;
 t(6)=30;
 %t(7,8)=[0.3 0.5];% fraction of ToF at which the DSM occurs
 t(7)=0.3;
 t(8)=0.5;
 t(9) =100;% perigee fly-by radius for planets P2..Pn-1, non-dimensional wrt planetary radii
 t(10) =1;% rotation gamma of the bplane-component of the swingby outgoing velocity (v_rel_out) [take n_r=cross(v_rel_in,v_planet_helio) if you rotate n_r by +gamma around v_rel_in you obtain the projection of v_rel_out on the b-plane] Vector (Vout) around the axis of the incoming swingby velocity vector (Vin)
 
MGADSMproblem.sequence=[3 3 4];%:  planet sequence. Example [3 3 4]= [E E M]. A     negative sign represents a retrograde（后向） orbit after the fly-by
MGADSMproblem.objective.type='orbit insertion';%: type of objective function (e.g. 'orbit insertion','rndv','gtoc1')
MGADSMproblem.objective.rp=100;%: pericentre radius of the target orbit if type = 'orbit insertion'
MGADSMproblem.objective.e=0;%:  eccentricity of the target orbit if type = 'orbit insertion'
MGADSMproblem.bounds=[10 1];%: decision vector upper and lower bounds for the optimiser
MGADSMproblem.yplot=1;%:     1-> plot trajectory, 0-> don't

[J,DVvec,DVarr] = GTOC6_mga_dsm(t,MGADSMproblem)




 