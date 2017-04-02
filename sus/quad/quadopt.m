function[pend] = quadopt


% coordinates x = longitudinal = u_LIGO roll about this axis

%********************************************************************
% MASS (N) is REPRESENTED by a rectangular BLOCK
% in reality it will be larger and less dense.
global pend
 pend.g		= 9.81;
 %pend.Inx	=  pend.mn*( pend.ny^2+ pend.nz^2)/12;	%moment of inertia (transverse roll)
% MASS (1) is REPRESENTED by a rectangular BLOCK
% in reality it will be larger and less dense.

 %pend.ux	= 0.13;			%dimensions MASS (1) (square)
 %pend.m1	= pend.den1* pend.uy* pend.uz* pend.ux	;	%MASS (1) 
%********************************************************************
 %pend.ix	= 0.13;			%dimension of MASS (2) (cylinder)
 pend.I2x	= pend.m2*(pend.ir^2/2);				%moment of inertia (transverse roll)
%*********************************************************************
 %pend.tx		= 0.13;			%dimensions of MASS (3) (cylinder)
 pend.I3x	= pend.m3*(pend.tr^2/2);				%moment of inertia (transverse roll)
%********************************************************************
%pend.ln = 0.54;    	%wire length 1: short version
%******************************************************************************************
 pend.nwn	= 2; %number of wires (= number of cantilevers if fitted) per stage (2 or 4)
%***********************************************************************************
%************************************************************************************
 %pend.Yn  	= 1.65e11;	%Youngs Modulus of wire	(N)	(s/steel 302)
 %blade design
 mntb = (pend.mn+pend.m1+pend.m2+pend.m3)/2;%total per blade
 %blade design
 m1tb = (pend.m1+pend.m2+pend.m3)/2;%total per blade
%blade design
 m2tb = (pend.m2+pend.m3)/2;%total per blade
%****************************************************************************************
%dees = 0.001;
pend.dm		= 0.001;	%height of wire break-off above c.of m. mass (N)
%additional information for ribbon breakoffs
%these are needed in translation/roll mode calculations
pend.twistlength = 0.00; %length of twist section in ribbon
%******************************************************************************************
%X direction separation 
pend.sn	= 0.00;		%1/2 separation of wires (N)   
%******************************************************************************************
%Y direction separation
pend.nn0 = 0.25;   	%1/2 separation of wires (N) at suspension point
%now we can work out the true lengths com to com etc.
pend.tln = sqrt(pend.ln^2 - (pend.nn0-pend.nn1)^2) + pend.dm;
%***********************************************************************************
% represents small loss
pend.bd = 0.00; % makes phases of open loop plots look nicer
%*********************************************************************************** 
%mean distances of actuators from axis of rotation
%pend.lever_pitch	= 0.06;