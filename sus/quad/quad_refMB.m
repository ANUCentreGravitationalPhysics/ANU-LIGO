% baseline referemce script for LIGO II MC/Recycling mirror suspension

% Uses code from Ken Strain 
% Calum I. Torrie November 2001
% this allows for only the parameters and the mode frequencies to be listed for a triple pendulum
% QUAD pendlulum parameters and frequencies


clear all
global pend

%ssmake4pv2e;%outdated, no longer supplied

%***********************************
%ssmake4pv2eMB; %updated files with errors derived from MATHMATICA, Mark Barton
%***********************************

ssmake4pv2eMB2; %further updated by Mark NArton to more correctly model blades+wires Feb05

   longpitch = ksort(damp(lpe)./2./pi);
   pend.longpitch1 = longpitch(1:4);
   pend.longpitch2 = longpitch(5:8);
   pend.yaw       = ksort(damp(ye)./2./pi);
   transroll = ksort(damp(trer)./2./pi);
   pend.transroll1 = transroll(1:4);
   pend.transroll2 = transroll(5:8);
   pend.vertical  = ksort(damp(ve)./2./pi)
   

