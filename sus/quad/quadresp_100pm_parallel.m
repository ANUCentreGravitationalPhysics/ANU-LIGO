% QUADRESP makes some plots of the quad pendulum response

% Assumes 'generate_simulink.m' has been run to generate the state-space
% model, 'quad_sys', of the quad pendulum

clear all;

% flag to save the figures
% 1 = yes, 0 = not
save_figure = 0;

f_low = 0.01;
f_high = 100;

f = logspace(log10(f_low),log10(f_high),1000);
w = 2*pi.*f;

% Angular radiation pressure torque coefficients
k_major = 0;
k_minor = 0;

% damper = 1; % ECD
% damper = 2; % GEO Damping
% damper = 3; % Damping with fancy LPF
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Transfer Function Suspension Point/GND -> TM
[mag(1,:), phs(1,:)] = bode(quad_sys(1,1),w);
H_GND_TM = mag .* exp(-i.*phs*pi/180);

clear mag phs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Feedback Servo Response, Open-loop

% PM Servo
fc_pm = 30;          % cut-off frequency, [Hz]
[zb, pb, kb] = besself(4,fc_pm*2*pi);
filter = zpk(zb, pb, kb);
%figure(99);
%bode(servo_filter, w);
[mag(1,:), phs(1,:)] = bode(filter,w);
cfilter = mag .* exp(i.*phs.*pi/180);

zz = [];
pp_pm = [1e-2];
kk = 20;
servo = zpk(-2*pi*zz, -2*pi*pp_pm, kk);
[mag(1,:), phs(1,:)] = bode(servo,w);
S_PM = mag .* exp(i.*phs.*pi/180) .* cfilter ./ H_PM_TM;

dummy = find(f > 10);
index_10 = dummy(1);
disp(['PM Servo gain at ',num2str(f(index_10),3),' Hz: ', num2str(20*log10(mag(index_10)),3), ' dB.']);
gainmargin = 20*log10(abs(mag(index_10)) - abs(H_PM_TM(index_10)));
clear mag phs;

% UIM Servo
fc_uim = 23;          % cut-off frequency, [Hz]
[zb, pb, kb] = besself(6,fc_uim*2*pi);
filter = zpk(zb, pb, kb);
%figure(99);
%bode(servo_filter, w);
[mag(1,:), phs(1,:)] = bode(filter,w);
cfilter = mag .* exp(i.*phs.*pi/180);

zz = [];
pp_uim = [1e-2];
kk = 0.9;
servo = zpk(-2*pi*zz, -2*pi*pp_uim, kk);
[mag(1,:), phs(1,:)] = bode(servo,w);
S_UIM = mag .* exp(i.*phs.*pi/180) .* cfilter ./ H_UIM_TM;

disp(['UIM Servo gain at ',num2str(f(index_10),3),' Hz: ', num2str(20*log10(mag(index_10)),3), ' dB.']);
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
X_PRN_noise = 100e-12 .* ones(size(f));  % flat displacement response at 100pm

% Test Mass displacement with SEI platform with a 
% 'standard' seismic level
% X_TM_sei = X_sei_std .* H_GND_TM;
X_GND_TM = X_sei_std .* H_GND_TM;

% Test Mass displacement with SEI platform and with a 
% 'noisy' seismic level
X_GND_TM_noisy = X_sei_noisy .* H_GND_TM;

% Test Mass displacement with SEI platform with a 'standard'
% seismic backgound and an SPI at 100pm at 1 Hz.
X_TM_goal = sqrt(1./((1./X_GND_TM).^2 + (1./X_PRN_noise).^2));


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close-Loop Displacements
G_PM = S_PM .* H_PM_TM;
G_UIM = S_UIM .* H_UIM_TM;

X_RO_noise = X_PRN_noise + X_GND_TM ./ (1 - G_PM - G_UIM);

% PM Feedback only
%X_TM = (X_GND_TM - X_PRN_noise .* G_PM) ./ (1 + G_PM);
%X_TM_rms = rms(f, abs(X_TM));

%X_PM_TM = (X_PRN_noise .* G_PM + X_GND_TM .* G_PM) ./ (1 + G_PM);
%X_PM_TM_rms = rms(f, abs(X_PM_TM));

% UIM + PM Feedback
X_TM = (X_GND_TM - X_PRN_noise .* G_UIM - X_PRN_noise .* G_PM) ./ ...
            (1 + G_UIM + G_PM);
X_TM_rms = rms(f, abs(X_TM));

X_UIM_TM = X_RO_noise .* G_UIM;
X_UIM_TM_rms = rms(f, abs(X_UIM_TM));

X_PM_TM = X_RO_noise .* G_PM;
X_PM_TM_rms = rms(f, abs(X_PM_TM));

%%%%%%%%%%%%%%%%%%%%%%%%%
% Feedback Force To Stage
FF_UIM = abs(X_RO_noise .* S_UIM);
FF_UIM_rms = rms(f, FF_UIM);

FF_PM = abs(X_RO_noise .* S_PM);
FF_PM_rms = rms(f, FF_PM);

FF_TM = abs(X_TM ./ H_TM_TM);
FF_TM_rms = rms(f, FF_TM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close-Loop Servo Response
HCL_SERVO_PM = X_PM_TM ./ X_RO_noise;
HCL_SERVO_UIM = X_UIM_TM ./ X_RO_noise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OSEM LSC Feedback Current
PM_OSEM = 4 * 26e-3;   %  4x Initial LIGO OSEM, Force coefficient, [N/A]
UIM_OSEM = 3 * 2.05;   % 3x AdvLIGO OSEM, Force Coefficient, [N/A]

I_OSEM_PM_max = 150e-3;     % Max Acquisition Current
I_OSEM_UIM_max = 10e-3;     % Max Acquisition Current

F_OSEM_PM_max = PM_OSEM * I_OSEM_PM_max;    % Max Total OSEM Force
F_OSEM_UIM_max = UIM_OSEM * I_OSEM_UIM_max; % Max Total OSEM Force

% Current to the individual PM OSEMs
I_OSEM_PM = FF_PM / PM_OSEM;
I_OSEM_UIM = FF_UIM / UIM_OSEM;

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
    print -dpng 100pm/f100_Quad_Transfer_Function.png
end

figure(110)
subplot(2,1,1)
hndl = semilogx(f, 20*log10(abs(S_PM)), ...
                f, 20*log10(abs(S_UIM)));
set(hndl, 'LineWidth', 2);
title('Servo Response - Servo');
ylabel('Magnitude [dB]');
axis([1e-2 100 60 180]);
grid on;
legend(['S_{PM}: Poles@',num2str(pp_pm),'Hz, 1/H_{PM->TM}, 4^{th}-order Bessel at ',num2str(fc_pm),'Hz'], ...
        ['S_{UIM}: Poles@',num2str(pp_uim),'Hz, 1/H_{UIM->TM}, 6^{th}-order Bessel at ',num2str(fc_uim),'Hz'], ...
        'Location', 'NorthWest');
subplot(2,1,2)
hndl = semilogx(f, angle(S_PM).*180/pi, ...
                f, angle(S_UIM).*180/pi);
set(hndl, 'LineWidth', 2);
ylabel('Phase [deg]');
xlabel('Frequency [Hz]');
axis([1e-2 100 -180 180]);
grid on;

if save_figure == 1
    print -dpng 100pm/f110_Servo_Response.png
end

figure(120)
hndl = loglog(f, abs(X_GND_TM), '-b', ...
              f, abs(X_TM_goal), '-g', ...
              f, abs(X_TM), '-r', ...
              f, abs(X_UIM_TM), 'm', ...
              f, abs(X_PM_TM), 'c', ...
              f, X_TM_rms, '--r');
axis([1e-2 100 1e-18 1e-5]);
grid on;
set(hndl, 'LineWidth', 2);
title('Displacement Noise - X');
ylabel('Test Mass Displacement Noise [m/rtHz]');
xlabel('Frequency [Hz]');
legend('TM without SPI', ...
       'TM-goal with PRN (100pm/rtHz)', ...
       'TM feedback (UIM+PM)', ...
       'TM with UIM', ...
       'TM with UIM + PM', ...
       'Location', 'SouthWest');

if save_figure == 1
    print -dpng 100pm/f120_TM_Displacement_Noise.png
end

figure(140)
hndl = loglog(f, I_OSEM_PM, '-b', ...
      f, I_OSEM_UIM, '-r', ...
      f, rms(f, I_OSEM_UIM), '--r', ...
      f, rms(f, I_OSEM_PM), '--b');
set(hndl, 'LineWidth', 2);
grid on;
axis([1e-2 1e2 1e-8 1e0]);
title('Single OSEM Feedback Current');
ylabel('Current [A/rtHz]');
xlabel('Frequency [Hz]');
legend('PM OSEM', ...
       'UIM OSEM', ...
       'Location', 'SouthWest');

if save_figure == 1
    print -dpng 100pm/f140_Single_OSEM_fb_Current.png
end

figure(150)
hndl = loglog(f, FF_PM, '-b', ...
              f, FF_TM, '-g', ...
              f, FF_UIM, '-r', ...
              f, FF_PM_rms, '--b', ...
              f, FF_TM_rms, '--g', ...
              f, FF_UIM_rms, '--r');
set(hndl, 'LineWidth', 2);
grid on;
axis([1e-2 1e2 1e-10 1e-2]);
title('Total Force To Stage');
ylabel('Force [N/rtHz]');
xlabel('Frequency [Hz]');
legend('F_{PM}', ...
       'F_{TM}', ...
       'F_{UIM}', ...
       'Location', 'SouthWest');

if save_figure == 1
    print -dpng 100pm/f150_Total_Force2Stage.png
end


figure(125)
subplot(2,1,1)
hndl = semilogx(f, 20*log10(abs(HCL_SERVO_PM)), ...
                f, 20*log10(abs(HCL_SERVO_UIM)));
set(hndl, 'LineWidth', 2);
title('Close-Loop Feedback Response - H = X_{RO,noise} / X_{..->TM}');
ylabel('Magnitude [dB]');
axis([1e-2 100 -80 100]);
grid on;
legend('H_{PM}', ...
       'H_{UIM}', ...
       'Location', 'SouthWest');
subplot(2,1,2)
hndl = semilogx(f, angle(HCL_SERVO_PM).*180/pi, ...
                f, angle(HCL_SERVO_UIM).*180/pi);
set(hndl, 'LineWidth', 2);
ylabel('Phase [deg]');
xlabel('Frequency [Hz]');
axis([1e-2 100 -180 180]);
grid on;

if save_figure == 1
    print -dpng 100pm/f125_Servo_Response_CL.png
end
