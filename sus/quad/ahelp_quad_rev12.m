% FILENAME:    ahelp_quad_rev12.m 
% DESCRIPTION: QUAD PENDULUM

% VERSION 1 Calum Torrie Jan 01
% VERSION 2 Phil Willems RM20020913
% VERSION 3 by CIT on 29th January 2003  
% VERSION 4 by CIT on 4th February 2003 
% VERSION % by CIT on 10th FEB 2003
% VERSION 9 by NAR on March 5th 2003
% VERSION 10 by CIT May 7th 2003
% VERSION 11 by CIT for a QUAD
% The generic names for the parameter file and pendulum file have been chosen so as no 
% changes are required in other files that call this set if the set of files is used
% for a new suspension.

% TEXT FILE EXPLAINING THE VARIOUS FILES ASSOCIATED WITH THE 
% ADVANCED LIGO QUADRUPLE PENDULUM MODEL

1) % TO OBTAIN PARAMETERS AND FREQUENCIES ONLY: -

 a) View and edit parameters of quadruple pendulum in 'quadopt.m'.
 
 b) Run 'quad_ref.m' for a list of parameters and frequencies. 'quad_ref.m' calls the spring
     matrix components from 'ssmake4pv2eMB.m' and calculates the blade parameters using 'opt.m'. 
     Opt.m should only be used as a starting guide to a blade design.
 
        
2) % TO OBTAIN CLOSED LOOP TRANSFER FUNCTIONS, IMPULSE RESPONSE etc..:-

 a) Again view and edit parameters in 'quadopt.m'.
 
 b) Open 'generate_simulink.m'. Here there a choice for either active damping or eddy current
     damping to be implemented. Set damper = 1 or 2  as required. There is a file called 
     'localdamp.m' which implements the appropriate damping calling 'local4ss.m' if required 
     (FOR ACTIVE DAMPING).
 
 c) Run 'generate_simulink.m'. This in turns calls the spring matrix components in 'ssmake3pv2n.m' 
     and the appropriate damping set in 2b). 
     'Ssmake3pv2n.m' calls 'Triplep.m' which in turn calls 'opt.m'. 
 
 d) 'Localdamp.m' calls 'highpass.m', 'lowpass.m' and 'transdif.m'.
    
    
3) % USING THE MODEL

 a) In PENDN under "tools" RUN "linear analysis".
 
 b) In LTI Viewer:PENDN under "Simulink" RUN "Get Linearised model".
 
 c) A reference of how the numbers are obtained for the gain triangle in PEND can be found at 
     the end of 'triplep.m'.
 
 d) A list of useful defaults for the user preferences can be found at the end of 'triplep.m'.


4) % USEFUL NOTES and FILES

  a) In 'quadopt.m' it is possible to input parameters that match either a "noise prototype" 
     or a "controls prototype" by replacing the appropriate section relating to the relevant
     mass and wire properties with a new set of numbers.
     A common practise used would be to have for example 2 sets of numbers for the upper mass
     with one set commented out using %%.

  b) AppendixB_quad.doc - layout explaining all of the key input parameters used in 'quadopt.m'

                    
5) % FILES NO LONGER CALLED or REQUIRED (These can be found in /outdated)

a) 'quadoptold.m' is an older input file.
b) 'qp_optimise.m', 'locl4s.m', 'pende.mdl' and 'lp_1' are left over from alternative methods
     used to implement the local control.
c) 'ssmake4pv2e.m' was previously used for eddy current damping and also to measure frequencies.
        
!!NOT INLCUDED!!   
   
6) % EXTRA FILES, USED ONLY AS ONE PROGRESSES TO THE FINE DETAIL OF THE DESIGN
%
% a)  'mofi2.m' Is a script that calculates the mass and moment of inertia of a "T" shaped upper mass.
% b)  'blades_etc.xls' can be used to calculate the following: -
%  (i)   Upper and lower cantilever blade parameters. 
%         (The calculation of the blade should again be treated with care, as this has to be used 
%          in conjunction with an ANSYS model of a blade.)
%  (ii)  Mass and moments of inertia of an intermediate mass with holes and flats. 
%  (iii) Mass and moments of inertia of a test mass with flats.
%  (iv)  Radius of wire to be used.       
