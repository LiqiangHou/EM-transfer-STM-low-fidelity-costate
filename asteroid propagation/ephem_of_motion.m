function [ephem_tf] = ephem_of_motion(ephem,tofdays)
% ephem:         [Index , ID, a_km,e,i_rad,o_rad,w_rad,M_rad,rp_au,ra_au,MJD_0 ]
% ephem_tf:      [Index , ID, a_km,e,i_rad,o_rad,w_rad,M_rad,rp_au,ra_au,MJD_f ]
% pos: km, 
% vel:km/s, 
% tof:  secs
global muSUN auKM daysecs
global min_tofdays

[n_ephem,n_dim] = size(ephem);


% tofsecs
if(tofdays < 0)
    tofdays = min_tofdays;
end
mjd_f   = ephem(1,11) + tofdays;
tofsecs = tofdays*daysecs;

for i=1:n_ephem

    a_km      = ephem(i,3);
    M_rad     = ephem(i,8);

    n_rad_sec = sqrt(muSUN / a_km^3);
    M_tf      = (M_rad + n_rad_sec*tofsecs);


    ephem_tf(i,:)   = [ephem(i,1:7), M_tf ,ephem(i,9:10), mjd_f];
end


return
end










