%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P. Fritschel 7 Aug 2001                                               %
% Script to generate model of pendulum and local damping function, and  %
% open the simulink model. Local damping functions are defined in       %
% 'localdamp.m' (can represent active or eddy current damping)          %
% added damping options K Strain 28 Aug 2001                            %
% new matrix elements for better blade modeling, M. Barton, 24 Feb 2005 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
global pend
pend.title = 'Pendulum parameters and derived properties';

%ssmake4pv2e; % Old, with errors
%ssmake4pv2eMB1; % Errors corrected from Mathematica, Mark Barton
%above line can be safely uncommented to compare old and new matrices.

%***********************************
ssmake4pv2eMB2; % better blade modeling from Mathematica, Mark Barton
%***********************************

pend      % print pendulum parameters

%set damper to correct value for desired damping option

%damper = 1;  %     eddy current
damper = 2;  %    adapted geo active

localdamp; 

simulink
open pende.mdl