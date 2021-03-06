function [ss_sys] = linquad2(dobie)

% LINQUAD makes a linearized state-space model of the quad pendulum
%
% 
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P. Fritschel 7 Aug 2001                                               %
% Script to generate model of pendulum and local damping function, and  %
% open the simulink model. Local damping functions are defined in       %
% 'localdamp.m' (can represent active or eddy current damping)          %
% added damping options K Strain 28 Aug 2001                            %
% new matrix elements for better blade modeling, M. Barton, 24 Feb 2005 %
%
% Changed structure to accomodate new .mdl file and new damping and
% noise source hooks.  Rana, Jul 2005
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
%clear all
global pend
pend.title = 'Pendulum parameters and derived properties';

%ssmake4pv2e; % Old, with errors
%ssmake4pv2eMB1; % Errors corrected from Mathematica, Mark Barton
%***********************************
ssmake4pv2eMB2; % better blade modeling from MATHMATICA, Mark Barton
%***********************************

%pend      % print pendulum parameters

%set damper to correct value for desired damping option

%damper = 1;  %     eddy current
%damper = 2;  %    adapted geo active
%damper = 3;  %   fancy low pass filtering
%damper = 4;  %   no damping

damper = dobie;

localdamp; 

% simulink
% open pende.mdl

% This makes a linearized state-space model of the
% quad pendulum Simulink model
[AAAA,BBBB,CCCC,DDDD] = linmod2('pende');
ss_sys = ss(AAAA,BBBB,CCCC,DDDD);
damper


