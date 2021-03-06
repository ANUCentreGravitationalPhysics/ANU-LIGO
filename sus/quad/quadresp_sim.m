% QUADRESP makes some plots of the quad pendulum response

% Assumes 'generate_simulink.m' has been run to generate the state-space
% model, 'quadResp', of the quad pendulum

clear all;

f_low = 0.01;
f_high = 100;

f = logspace(log10(f_low),log10(f_high),1000);
w = 2*pi.*f;

% Angular radiation pressure torque coefficients
k_major = 0;
k_minor = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Feedback Servo Response, Open-loop

fc_bessel = 25;         % filter cut-off frequency [Hz]
[zCutOff, pCutOff, kCutOff] = besself(4,fc_bessel*2*pi);

CutOffFilter = zpk(zCutOff, pCutOff, kCutOff);
%figure(99);
%bode(servo_filter, w);
[mag(1,:), phs(1,:)] = bode(CutOffFilter,w);
cfilter = mag .* exp(i.*phs.*pi/180);
clear mag phs;

% PM Servo
zPM = [0.1 0.4 1.8 1.8 zCutOff/2/pi];
pPM = [1e-3 pCutOff'/2/pi];
kPM = 700 * kCutOff;

% Set the Simulink Blocks
zpkPM = zpk(-2*pi*zPM, -2*pi*pPM, kPM);
% bode(servo_pm);

[mag(1,:), phs(1,:)] = bode(zpkPM,w);
S_PM = mag .* exp(i.*phs.*pi/180); 

% % UIM Servo
% zz = [];
% pp = [1e-2];
% kk = 200;
% servo = zpk(-2*pi*zz, -2*pi*pp, kk);
% [mag(1,:), phs(1,:)] = bode(servo,w);
% S_UIM = mag .* exp(i.*phs.*pi/180) .* cfilter ./ H_UIM_TM;

% dummy = find(f > 10);
% index_10 = dummy(1);
% disp(['Servo gain at ',num2str(f(index_10),3),' Hz: ', num2str(20*log10(mag(index_10)),3), ' dB.']);
% gainmargin = 20*log10(abs(mag(index_10)) - abs(H_PM_TM(index_10)))

clear mag phs;

% Set the Simulink Blocks
ServoPM = zpkPM;
ServoUIM = 0;

% Switches to switch between Open and Close loop response
% switch > 0 - Close Loop Response
% switch < 0 - Open Loop Response
PMswitch = -1;
UIMswitch = -1;

% damper = 1; % ECD
% damper = 2; % GEO Damping
% damper = 3; % Damping with fancy LPF
% damper = 4; % no damping
damper = 3;
linquadLAI;
quadOL = quad_sys;
%save quadsys quadResp

% Close Loop Respnse
PMswitch = 1;
UIMswitch = 1;
linquadLAI;
quadCL = quad_sys;

% quadResp = quadCL;
quadResp = quadCL;

%%%%%%%%%%%%%%%%%%%%
%  Transfer Function
% pende model
% Force input               Displ. output
% 1-GND                     1- X-TM
% 18-TOP                    7- X-PM
% 10-UIM                    8- X-UIM
% 9-PM                      9- PM Feedback Force
% 8-TM                      10- UIM Feedback Force
% 21- PRN Noise Limit       11- PM Servo Input
% 22- OL PM Servo Input
% 23- OL UIM Servo Input
%
% syntax: quadResp(output, input)

[mag(1,:), phs(1,:)] = bode(quadResp(1,8),w);
H_TM_TM = mag .* exp(i.*phs.*pi/180);

[mag(1,:), phs(1,:)] = bode(quadResp(1,9),w);
H_PM_TM = mag .* exp(i.*phs.*pi/180);

[mag(1,:), phs(1,:)] = bode(quadResp(7,10),w);
H_UIM_PM = mag .* exp(i.*phs*pi/180);

[mag(1,:), phs(1,:)] = bode(quadResp(1,10),w);
H_UIM_TM = mag .* exp(i.*phs*pi/180);

[mag(1,:), phs(1,:)] = bode(quadResp(8,18),w);
H_TOP_UIM = mag .* exp(i.*phs*pi/180);

clear mag phs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Transfer Function Suspension Point/GND -> TM
[mag(1,:), phs(1,:)] = bode(quadResp(1,1),w);
H_GND_TM = mag .* exp(i.*phs*pi/180);

[mag(1,:), phs(1,:)] = bode(quadResp(7,1),w);
H_GND_PM = mag .* exp(i.*phs*pi/180);

clear mag phs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Transfer Function PM Servo -> TM
%  This is also the openloop TF
[mag(1,:), phs(1,:)] = bode(quadResp(1,21),w);
H_PRN_TM = mag .* exp(i.*phs*pi/180);

[mag(1,:), phs(1,:)] = bode(quadResp(11,21),w);
H_PRN_PMservo = mag .* exp(i.*phs*pi/180);

clear mag phs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Transfer Function TM -> PM Force Feedback
%[mag(1,:), phs(1,:)] = bode(quadResp(9,8),w);
%HF_TM_PM = mag .* exp(i.*phs*pi/180);
%
%clear mag phs;
[mag(1,:), phs(1,:)] = bode(quadResp(9,21),w);
H_PRN_Fpm = mag .* exp(i.*phs*pi/180);

[mag(1,:), phs(1,:)] = bode(quadResp(9,1),w);
H_GND_Fpm = mag .* exp(i.*phs*pi/180);

clear mag phs;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Open-loop Transfer Function
% G_PM = S_PM .* H_PM_TM;
H_GND_PRN = H_GND_TM ./ H_PRN_TM;
H_OL = H_GND_PRN;

H_PRN_PM = H_PRN_TM ./ H_PM_TM;

H_TM_Fpm = H_GND_Fpm ./ H_GND_TM;
H_CL = H_TM_Fpm;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Generate the Seismic Platform Displacement Requirement
seismic_platform_displ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPI/LAI PRN displacement noise,
% flat displacement response at 100pm
X_PRN_noise = 100e-12 .* ones(size(f));  

% Test Mass displacement with SEI platform with a 
% 'standard' seismic level
X_GND_TM = X_sei_std .* abs(H_GND_TM);
X_GND_TM_rms = rms(f, abs(X_GND_TM));

% Test Mass displacement with SEI platform and with a 
% 'noisy' seismic level
% X_GND_TM_noisy = X_sei_noisy .* H_GND_TM;

% Test Mass displacement with SEI platform with a 'standard'
% seismic backgound and an SPI at 100pm at 1 Hz.
X_TM_goal = sqrt(1./((1./X_GND_TM).^2 + (1./X_PRN_noise).^2));


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close-Loop Displacements
X_PRN_TM = abs(H_PRN_TM) .* X_GND_TM;
% X_TM_PRN_rms = rms(f, abs(X_TM_PRN));
% 
% X_PM_TM = (X_PRN_noise .* G_PM + X_GND_TM .* G_PM) ./ (1 + G_PM);
% X_PM_TM_rms = rms(f, abs(X_PM_TM));

%%%%%%%%%%%%%%%%%%%%%%%%%
% Feedback Force To Stage
% Term1 = X_PRN_noise .* S_PM;
% Term2 = X_GND_TM .* S_PM;
%Term2 = zeros(size(f));

% FF_PM = (Term1 + Term2) ./ (1 + G_PM);
% FF_PM = abs(FF_PM);
% FF_PM_rms = rms(f, FF_PM);
% 
% FF_TM = abs(X_TM_PRN ./ H_TM_TM);
% FF_TM_rms = abs(X_TM_PRN_rms ./ H_TM_TM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close-Loop Servo Response
% X_RO_noise = X_PRN_noise + X_TM_sei ./ (1 + G_PM);
% 
% H_SERVO_OL1 = G_PM;
% 
% H_SERVO_OL2 = X_PM_TM ./ X_RO_noise;
% 
% H_SERVO_CL = X_TM_PRN ./ X_GND_TM;
% 
% H_SUPP = X_RO_noise ./ X_TM_PRN;
% 
% H_RESID = X_PM_TM ./ X_TM_PRN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OSEM LSC Feedback Current
PM_OSEM = 4 * 26e-3;   %  4x Initial LIGO OSEM, Force coefficient, [N/A]
UIM_OSEM = 3 * 2.05;   % 3x AdvLIGO OSEM, Force Coefficient, [N/A]

I_OSEM_PM_max = 400e-3;     % Max Acquisition Current
I_OSEM_UIM_max = 10e-3;     % Max Acquisition Current

F_OSEM_PM_max = PM_OSEM * I_OSEM_PM_max;    % Max Total OSEM Force
F_OSEM_UIM_max = UIM_OSEM * I_OSEM_UIM_max; % Max Total OSEM Force

% Current to the individual PM OSEMs
I_OSEM_PM = abs(H_TM_Fpm) / PM_OSEM;
I_OSEM_PM_rms = rms(f, I_OSEM_PM);

quadresp_plot;
