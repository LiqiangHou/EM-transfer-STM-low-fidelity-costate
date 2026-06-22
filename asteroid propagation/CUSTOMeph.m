
%--------------------------------------------------------------------------
function [r,v]=CUSTOMeph(jd,epoch,keplerian,flag)
%
%Returns the position and the velocity of an object having keplerian
%parameters epoch,a,e,i,W,w,M
%
%Usage:     [r,v]=CUSTOMeph(jd,name,list,data,flag)
%
%Inputs:    jd: julian date
%           epoch: mjd when the object was observed (referred to M)
%           keplerian: vector containing the keplerian orbital parameters
%
%Output:    r = object position with respect to the Sun (km if flag=1, AU otherwise)
%           v = object velocity ( km/s if flag=1, AU/days otherways )
%
%Revisions :    Function added 04/07

global AU mu

muSUN = mu(11);
     a=keplerian(1)*AU; %in km
     e=keplerian(2);
     i=keplerian(3); 
     W=keplerian(4);
     w=keplerian(5);
     M=keplerian(6);
     jdepoch=mjd2jed(epoch);
     DT=(jd-jdepoch)*60*60*24;
     n=sqrt(muSUN/a^3);
     M=M/180*pi;
     M=M+n*DT;
     M=mod(M,2*pi);
     E=M2E(M,e);
     [r,v]=par2IC([a,e,i/180*pi,W/180*pi,w/180*pi,E],muSUN);
     if flag~=1
         r=r/AU;
         v=v*86400/AU;
     end
 

