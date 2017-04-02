% QUADRESP makes some plots of the quad pendulum response

% Assumes 'generate_simulink.m' has been run to generate the state-space
% model, 'quad_sys', of the quad pendulum

clear all;

% flag to save the figures
% 1 = yes, 0 = not
save_figure = 1;

f_low = 0.01;
f_high = 100;

f = logspace(log10(f_low),log10(f_high),1000);
w = 2*pi.*f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Feedback Servo Response, Open-loop
fc_bessel = 25;         % filter cut-off frequency [Hz]
[zb, pb, kb] = besself(4,fc_bessel*2*pi);
servo_filter = zpk(zb, pb, kb);
%figure(99);
%bode(servo_filter, w);
[mag(1,:), phs(1,:)] = bode(servo_filter,w);
cfilter = mag .* exp(i.*phs.*pi/180);

% PM Servo
zz = 1.5*[1 1 1 1];
pp = [1e-2];
kk = 30;
servo = zpk(-2*pi*zz, -2*pi*pp, kk);
[mag(1,:), phs(1,:)] = bode(servo,w);
% S_PM = mag .* exp(i.*phs.*pi/180) .* cfilter ./ H_PM_TM; % 1/H_PM_TM as
% servo filter!
S_PM = mag .* exp(i.*phs.*pi/180) .* cfilter;

% UIM Servo
zz = 1*[1 1 1 1];
pp = [1e-2];
kk = 200;
servo = zpk(-2*pi*zz, -2*pi*pp, kk);
[mag(1,:), phs(1,:)] = bode(servo,w);
S_UIM = mag .* exp(i.*phs.*pi/180) .* cfilter;

% dummy = find(f > 10);
% index_10 = dummy(1);
% disp(['Servo gain at ',num2str(f(index_10),3),' Hz: ', num2str(20*log10(mag(index_10)),3), ' dB.']);
% gainmargin = 20*log10(abs(mag(index_10)) - abs(H_PM_TM(index_10)))

clear mag phs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quad Simulink Model, with SPI/LAI Feedback

% Angular radiation pressure torque coefficients
k_major = 0;
k_minor = 0;

% LAI Feedback servo
ServoPM = S_PM;
ServoUIM = S_UIM;

% Quad Simulink Model
% Quad local dmping model
% damper = 1; % ECD
% damper = 2; % GEO Damping
% damper = 3; % Damping with fancy resoance and LPF
% damper = 4; % no damping
damper = 3;
linquad;
quad_sys = ss_sys;
%save quadsys quad_sys

%%%%%%%%%%%%%%%%%%%%
%  Transfer Function
[mag(1,:), phs(1,:)] = bode(quad_sys(1,8),w);
H_TM_TM = mag .* exp(-i.*phs.*pi/180);

[mag(1,:), phs(1,:)] = bode(quad_sys(1,9),w);
H_PM_TM = mag .* exp(-i.*phs.*pi/180);

[mag(1,:), phs(1,:)] = bode(quad_sys(7,10),w);
H_UIM_PM = mag .* exp(-i.*phs*pi/180);

[mag(1,:), phs(1,:)] = bode(quad_sys(1,10),w);
H_UIM_TM = mag .* exp(-i.*phs*pi/180);

[mag(1,:), phs(1,:)] = bode(quad_sys(8,18),w);
H_TOP_UIM = mag .* exp(-i.*phs*pi/180);

clear mag phs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Transfer Function Suspension Point/GND -> TM
[mag(1,:), phs(1,:)] = bode(quad_sys(1,1),w);
H_GND_TM = mag .* exp(-i.*phs*pi/180);

clear mag phs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seismic Platform Displacement Noise
% This is a noisy estimate (lots of people, say in 10 years time)
N_F = [0.01 0.03 0.1 0.15 0.2 0.3 0.5 1 10 30 100];
N_X = [5e-6 1e-6 1e-6 2e-6 5e-7 2e-7 8e-10 3e-11 6e-13 6e-14 3e-15];
sei_noisy_0612 = 10.^(interp1(N_F,log10(N_X),f,'cubic',-14));
X_sei_noisy = sei_noisy_0612;

% Seismic Platform Displacement Noise
% The is the current requirement
SEI_F = [0.01 0.03 0.1 0.2 0.5 1 10 30 100];
SEI_X = [3e-6 1e-6 2e-7 2e-7 8e-10 1e-11 3e-13 3e-14 1e-15];
sei_model_0612 = 10.^(interp1(SEI_F,log10(SEI_X),f,'cubic',-14));
X_sei_std = sei_model_0612;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPI/LAI PRN displacement noise, 100pm at 1 Hz (1/f - response)
% X_spi = 100e-12 ./ f;
X_PRN_noise = 3e-12 .* ones(size(f));  % flat displacement response at 100pm

% Test Mass displacement with SEI platform with a 
% 'standard' seismic level
X_TM_sei = X_sei_std .* H_GND_TM;
X_GND_TM = X_TM_sei;
X_GND_TM_rms = rms(f, abs(X_GND_TM));

% Test Mass displacement with SEI platform and with a 
% 'noisy' seismic level
X_TM_noisy = X_sei_noisy .* H_GND_TM;

% Test Mass displacement with SEI platform with a 'standard'
% seismic backgound and an SPI at 100pm at 1 Hz.
X_TM_goal = sqrt(1./((1./X_TM_sei).^2 + (1./X_PRN_noise).^2));


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close-Loop Displacements
G_PM = S_PM .* H_PM_TM;
X_TM_PRN = (X_GND_TM - X_PRN_noise .* G_PM) ./ (1 + G_PM);
X_TM_PRN_rms = rms(f, abs(X_TM_PRN));

X_PM_TM = (X_PRN_noise .* G_PM + X_GND_TM .* G_PM) ./ (1 + G_PM);
X_PM_TM_rms = rms(f, abs(X_PM_TM));

%%%%%%%%%%%%%%%%%%%%%%%%%
% Feedback Force To Stage
Term1 = X_PRN_noise .* S_PM;
Term2 = X_GND_TM .* S_PM;
%Term2 = zeros(size(f));

FF_PM = (Term1 + Term2) ./ (1 + G_PM);
FF_PM = abs(FF_PM);
FF_PM_rms = rms(f, FF_PM);

FF_TM = abs(X_TM_PRN ./ H_TM_TM);
FF_TM_rms = abs(X_TM_PRN_rms ./ H_TM_TM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close-Loop Servo Response
X_RO_noise = X_PRN_noise + X_TM_sei ./ (1 - G_PM);

H_SERVO_CL = X_PM_TM ./ X_RO_noise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OSEM LSC Feedback Current
PM_OSEM = 4 * 26e-3;   %  4x Initial LIGO OSEM, Force coefficient, [N/A]
UIM_OSEM = 3 * 2.05;   % 3x AdvLIGO OSEM, Force Coefficient, [N/A]

I_OSEM_PM_max = 400e-3;     % Max Acquisition Current
I_OSEM_UIM_max = 10e-3;     % Max Acquisition Current

F_OSEM_PM_max = PM_OSEM * I_OSEM_PM_max;    % Max Total OSEM Force
F_OSEM_UIM_max = UIM_OSEM * I_OSEM_UIM_max; % Max Total OSEM Force

% Current to the individual PM OSEMs
I_OSEM_PM = FF_PM / PM_OSEM;
I_OSEM_PM_rms = rms(f, I_OSEM_PM);

%%%%%%%%%%
% Plotting
figure(100)
hndl = loglog(f, abs(H_TM_TM), ...
       f, abs(H_PM_TM), ...
       f, abs(H_UIM_TM));
set(hndl, 'LineWidth', 2);
grid on;
title('Quad Transfer Function - H = x/F')
xlabel('Frequency [Hz]');
ylabel('Magnitude [m/N]');
legend('TM -> TM', ...
       'PM -> TM', ...
       'UIM -> TM', ...
       'Location', 'NorthEast');
axis([0.01 100 1e-8 1e0])

if save_figure == 1
    print -dpng 2pm/f100_Quad_Transfer_Function.png
end

figure(110)
subplot(2,1,1)
hndl = semilogx(f, 20*log10(abs(S_PM)));
set(hndl, 'LineWidth', 2);
title('Servo Response - S_{PM}');
ylabel('Magnitude [dB]');
axis([1e-2 100 60 180]);
grid on;
legend(['Poles@',num2str(pp),'Hz, Zeros@', num2str(zz),'Hz, 4^{th}-order Bessel at ',num2str(fc_bessel),'Hz'], ...
       'Location', 'NorthWest');
subplot(2,1,2)
hndl = semilogx(f, angle(S_PM).*180/pi);
set(hndl, 'LineWidth', 2);
ylabel('Phase [deg]');
xlabel('Frequency [Hz]');
axis([1e-2 100 -180 180]);
grid on;

if save_figure == 1
    print -dpng 2pm/f110_Servo_Response.png
end

figure(120)
hndl = loglog(f, abs(X_GND_TM), ...
              f, abs(X_PRN_noise), ...
              f, abs(X_TM_PRN), ...
              f, X_TM_PRN_rms, '--r', ...
              f, X_GND_TM_rms, '--b');
axis([1e-2 100 1e-18 1e-5]);
grid on;
set(hndl, 'LineWidth', 2);
title('Displacement Noise - X');
ylabel('Test Mass Displacement Noise [m/rtHz]');
xlabel('Frequency [Hz]');
legend('TM without SPI', ...
       ['PRN noise level (',num2str(X_PRN_noise(1),3),' m/rtHz)'], ...
       'TM with PRN feedback', ...
       'Location', 'SouthWest');

if save_figure == 1
    print -dpng 2pm/f120_TM_Displacement_Noise.png
end

figure(140)
hndl = loglog(f, I_OSEM_PM, ...
      f, I_OSEM_PM_rms, '--b');
set(hndl, 'LineWidth', 2);
grid on;
axis([1e-2 1e2 1e-8 1e0]);
title('Single OSEM Feedback Current');
ylabel('Current [A/rtHz]');
xlabel('Frequency [Hz]');
legend(['PM OSEMs - ',num2str(I_OSEM_PM_rms(1)*1000,3),' mA_{rms}'], ...
       'Location', 'SouthWest');

if save_figure == 1
    print -dpng 2pm/f140_Single_OSEM_fb_Current.png
end

figure(150)
hndl = loglog(f, FF_PM, '-b', ...
              f, FF_TM, '-g', ...
              f, FF_PM_rms, '--b', ...
              f, FF_TM_rms, '--g');
set(hndl, 'LineWidth', 2);
grid on;
axis([1e-2 1e2 1e-10 1e-2]);
title('Total Force To Stage');
ylabel('Force [N/rtHz]');
xlabel('Frequency [Hz]');
legend('F_{PM}', ...
       ['F_{TM} - ', num2str(FF_TM_rms(1)*1e6,3),' uN_{rms}'], ...
       'Location', 'SouthWest');

if save_figure == 1
    print -dpng 2pm/f150_Total_Force2Stage.png
end


figure(125)
subplot(2,1,1)
hndl = semilogx(f, 20*log10(abs(H_SERVO_CL)));
set(hndl, 'LineWidth', 2);
title('Close-Loop Feedback Response - H = X_{RO,noise} / X_{PM->TM}');
ylabel('Magnitude [dB]');
axis([1e-2 100 -50 100]);
grid on;
%legend(['Poles@',num2str(pp),'Hz, Zeros@', num2str(zz),'Hz'], ...
%       'Location', 'NorthWest');
subplot(2,1,2)
hndl = semilogx(f, angle(H_SERVO_CL).*180/pi);
set(hndl, 'LineWidth', 2);
ylabel('Phase [deg]');
xlabel('Frequency [Hz]');
axis([1e-2 100 -180 180]);
grid on;

if save_figure == 1
    print -dpng 2pm/f125_Servo_Response_CL.png
end
