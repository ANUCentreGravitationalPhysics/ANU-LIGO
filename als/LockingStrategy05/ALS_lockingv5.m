% ALS Locking Strategy
%
% Daniel Sigg's idea using 4 VCO's
%
% Needs the ALS_freq3v3.mdl Simulink model
%
% BS - 10 May 2010
%
% 24 June 2010 - modified the titles to reflect the G = b/a, e.g. (in 18, out 22/out 20)
% 13 July 2010 - Added more input/output points to get Hz/Hz transfer
%                functions for VCO noise performance requirements
% 17 July 2010 - Removed the cutoff filters in the TM and PM feedback, 
%                In creased the Comm BW to ~10 kHz
%                Increased the Diff feedback to ~55 Hz by increasing the
%                cutoff filters to 250 Hz and 150 Hz
% v5 - 15 Nov 2010 - included signal injection into the X-ground motion
%                into the quad.
% v5 - 11 Dec 2010 - included signal injection into the PSL SHG system to
%                to monitor the crystal temperature drift induced phase
%                change of the GRN w.r.t. the PSL IR. -> figure 24.
% v6 - 4 July 2011 - added velocity to the TM displacement plot (figure 22)
%

clear all;
close all;

% Constants
cc = 299792458; % [m/s]

lambda_IR = 1064e-9;
lambda_GRN = 532e-9;

% f = logspace(-2, 2, 1e3);  % used with the SEI transfer function
 f = logspace(-1, 5, 1e3);

%% Engaging Servo's
FSSengaged = 1;
TTFSSengaged = 1;
PDHengaged = 1;
COMMengaged = 1;
DIFFengaged = 1;

%% Setting up the PSL section
%

% Ref Cav Transmission TF
RefCavFSR = cc / (2 * 0.2); % Ref Cav length 20cm?
RefCavFIN = 5000; % Ref Cav Finesse
RefCavPole = 2*pi * RefCavFSR / (2* RefCavFIN); % Ref Cav Pole frequency, 470kHz

% PSL FSS
PSLfrequencyActuator = 1;     % This is just a gain for the FSS feedback to the PSL (Temp, PZT and Pockell)
zzz = [RefCavPole/2/pi];
ppp = [1];
kkk = 10^( 17 /20);
FSS_Servo = zpk(-2*pi*zzz, -2*pi*ppp, -kkk);

% LowNoiseVCO -> FSS low noise VCO driver TF, [Hz/V]
zzz = [];
ppp = 2e6;                  % Range of the VCO
kkk = ppp / 20;             % VCO Full tuning range (2 MHz) / VCO Input voltage range (+/-20V)
LowNoiseVCO_psl = zpk(-2*pi*zzz, -2*pi*ppp, kkk);

% Noisy Ref Cav locking @ 1kHz rms below 10 Hz
zzz = [];
ppp = 10;                  % 
kkk = 6.2e2 * ppp;           % 
RefCavNoisy = zpk(-2*pi*zzz, -2*pi*ppp, kkk);
NoisyRefCav = mybodesys(RefCavNoisy, f);
% figure(1)
% loglog(f, abs(NoisyRefCav), f, rms(f', abs(NoisyRefCav)));
% axis([1e-1 1e4 1e-1 1e3]);

% AOM Driver at fiber launch
AOM_Driver = 1;

%% Setting up the X-End station
%

FiberPhaseNoise = 100;   % flat fiber induced phase noise, 100 Hz/rtHz
freqnoiseNPRO = abs(1e4 ./ (1 + i.*f/1)); % Freerunning NPRO, 100 Hz/rtHz at 100 Hz

% Arm Cav Transmission TF
L_arm = 3995;
ArmCavFSR = cc / (2 * L_arm); % Ref Cav length 20cm?
ArmCavFIN = 100; % Ref Cav Finesse
ArmCavPole = 2*pi* ArmCavFSR / (2* ArmCavFIN); % Arm Cav Pole frequency, ~1178 Hz

% LowNoiseVCO,EX -> TTFSS low noise VCO driver TF, used to demodulate the
% heterodyne signal in the end-station
zzz = [];
ppp = 2e6;
kkk = ppp / 20;             % VCO Full tuning range (2 MHz) / VCO Input voltage range (+/-20V)
LowNoiseVCO_local = zpk(-2*pi*zzz, -2*pi*ppp, kkk);

% End-Station laser feedback
uPZT_local = 3e6;   % PZT Volts to IR Frequency conversion, 3 MHz/V

% TTFSS Servo -> TTFSS Locking Servo, lock the laser to the heterodyne beatnote.
zzz = [];%[ 0 1e3];%[0 0 10 RefCavPole/2/pi];
ppp = [2e5];%[1 1 2e5];            % limited by the feedback to the laser PZT?
kkk = 10^(80/20);%10^(79/20);
TTFSS_Servo = zpk(-2*pi*zzz, -2*pi*ppp, kkk);

% PDH Servo -> PDH locking servo to lock the laser frequency to the arm cavity
% add notch at 1 Hz, Q=10, depth 30 dB (using foton:)

%- When DIFF and COMM are not engaged
% zzz = [1 300 0.003+i*0.999949 0.003-i*0.999949];    
% ppp = [1.001+i*0.994872 1.001-i*0.994872 1e5 1e5];
% kkk = 100000000000000;%6000000;

zzz = [200 ];    
ppp = [1];
kkk = 10^(62/20);
PDH_Servo = zpk(-2*pi*zzz, -2*pi*ppp, kkk);

%% Setting up the Vertex ALS Demodulation
%

% LowNoiseVCO,Comm -> Vertex ALS Common Mode low noise VCO driver TF
zzz = [];
ppp = 2e6;                  % Frequency range, Hz
kkk = ppp / 20;             % VCO Full tuning range (2 MHz) / VCO Input voltage range (+/-20V)
LowNoiseVCO_comm = zpk(-2*pi*zzz, -2*pi*ppp, kkk);

% LowNoiseVCO,Diff -> Vertex ALS Differential Mode low noise VCO driver TF
LowNoiseVCO_diff = zpk(-2*pi*zzz, -2*pi*ppp, kkk);

%% Common Mode Servo -> Common Mode locking servo to lock the PSL frequency
%% to the common mode arm cavity length fluctuations
zzz = [0];
ppp = [1];
%kkk = 10^( 126 /20);
kkk = 10^( 190 /20);
CommonMode_Servo = zpk(-2*pi*zzz, -2*pi*ppp, -kkk);

%% Differential Signal Feedback to the ETM Quads
nu_IR = cc / lambda_IR;
nu_GRN = cc / lambda_GRN;

% Differential Mode Servo
zzz = [0.5];
ppp = [1e6];
%kkk = 10^( -50 /20);
kkk = 10^( -55 /20);
DifferentialMode_Servo = zpk(-2*pi*zzz, -2*pi*ppp, -kkk);

%% Setting up the Quad Feedback and Control Block
global pend

% Angular radiation pressure torque coefficients
k_major = 0;
k_minor = 0;
k_ospring = 0;

gLP = 0;
ServoTM = 0;
ServoPM = 0;
ServoUIM = 0;
ServoTOP = 0;

damper = 1; % ECD
% damper = 2; % GEO Damping
% damper = 3; % Damping with fancy LPF
% damper = 4; % no damping

%***********************************
ssmake4pv2eMB2; % better blade modeling from MATHEMATICA, Mark Barton
%***********************************
localdamp;

if ~exist('k_ospring')
  k_ospring = 0;
else
  % MEVANS
  warning('NOT including optical spring.');
  k_ospring = 0;
end

% Run the Quad Servo Script
% PDHservo_All_2010_05_21_11_00_19 % original with 1pm noise floor
PDHservo_All_2010_07_19_15_51_34  % new 0.01pm noise floor

% Set the Signal Path Switches
gLP = 0;    % Keep the loop open to make it run within the overall simulation
gTM = 1;    % engage the TM feedback
gPM = 1;    % engage the PM feedback
gUIM = 0;   % engage the UIM feedback
gTOP = 0;   % engage the TOP feedback

%% Implementing the Simulink Model
%
%modelname = 'ALS_freq3v3';
%modelname = 'ALS_freq3v4';
modelname = 'ALS_freq3v5';
%
[AAA,BBB,CCC,DDD] = linmod2(modelname); % linearise the Simulink model
[rw,cl] = find(AAA == Inf); AAA(rw,cl) = 1e20;
[rw,cl] = find(AAA == -Inf); AAA(rw,cl) = -1e20;
[rw,cl] = find(BBB == Inf); BBB(rw,cl) = 1e20;
[rw,cl] = find(BBB == -Inf); BBB(rw,cl) = -1e20;
[rw,cl] = find(CCC == Inf); CCC(rw,cl) = 1e20;
[rw,cl] = find(CCC == -Inf); CCC(rw,cl) = -1e20;
[rw,cl] = find(DDD == Inf); DDD(rw,cl) = 1e20;
[rw,cl] = find(DDD == -Inf); DDD(rw,cl) = -1e20;

SYS = ss(AAA,BBB,CCC,DDD);              % ceates a state-space model of the Simulink model
                                        % TF = SYS(output, input)

%% Do some test plotting
% a = SYS(27,25);
% a = SYS(4,25);
% b = SYS(31,4);
% figure(99)
% mybodesys(a,f);
%title('VCO\_EX -to- f\_IR\_ex');


%% Print the Simulink model with all its colours, that works only in
%% Windows, and Mac OS X with matlab version R2010b!
% set_param(modelname, 'ShowPageBoundaries', 'on');    
% print(['-s' modelname], '-dpdf', [modelname '.pdf']); % print the simulink model with its colors...
    
%% Obtaining the transfer functions
save_figure = 0;    % controls is the figures are save as .pdf or not

save_figure_all = 0;
save_figure_dir = 'sim/';

if save_figure_all
    FSS=1;
    TTFSS=1;
    PDH=1;
    COMM=1;
    DIFF=1;
    save_figure = 1;
    save_figure_dir = 'sim/';
else
    FSS = 0;    % Plots the FSS loops of th PSL servo
    if FSSengaged
        FSS = 1;
        save_figure_dir = 'ttfss/';
    end
    TTFSS = 0;  % Plots the TTFSS loops in the end-station
    if TTFSSengaged
        TTFSS = 1;
        save_figure_dir = 'ttfss/';
    end
    PDH = 0;    % plots the PDH loops in the end-station
    if PDHengaged
        PDH = 1;
        save_figure_dir = 'pdh/';
    end
    COMM = 0;
    if COMMengaged
        COMM = 1;
        save_figure_dir = 'comm/';
    end
    DIFF = 0;
    if DIFFengaged
        DIFF = 1;
        save_figure_dir = 'diff/';
    end
end         % end save_figure

if FSS
%% FSS Feedback of the laser to the Reference Cavity
    %
%%  % PSL laser -to- FSS error signal
    hdl= figure(101)
    a = SYS(20,18);
    b = SYS(22,18);
    G = mybodesys(b/a,f);
    %
    subplot(211)
    semilogx(f,20*log10(abs(G)), 'LineWidth', 2)
    tt= title('PSL -to- FSS error signal TF (in 18, out 22/out 20)', 'FontSize', 16);
    ylabel('Mag [dB]', 'FontSize', 14);
    axis([min(f) max(f) floor(min(log10(abs(G))))*20 ceil(max(log10(abs(G))))*20])
    grid on
    subplot(212)
    semilogx(f,angle(G)*180/pi, 'LineWidth', 2)
    ylabel('Phase [deg]', 'FontSize', 14);
    xlabel('Freq [Hz]', 'FontSize', 14);
    axis([min(f) max(f) 90*floor(min(angle(G)*180/pi)/90) 90*ceil(max(angle(G)*180/pi)/90)])
%    set(gca,'YTick',[-1800:90:1800])
    grid on
    if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
        % print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), ' ',get(tt,'string'), '.pdf']);
    end


%%  % FSS Open Loop response
    hdl= figure(102)
    a = SYS(21,19);
    b = SYS(20,19);
    G = mybodesys(b/a,f);
    %
    subplot(211)
    semilogx(f,20*log10(abs(G)), 'LineWidth', 2)
    tt= title('FSS Controller TF (in 19, out 20/out 21)','FontSize', 16);
    ylabel('Mag [dB]', 'FontSize', 14);
    axis([min(f) max(f) floor(min(log10(abs(G))))*20 ceil(max(log10(abs(G))))*20])
    grid on
    subplot(212)
    semilogx(f,angle(G)*180/pi, 'LineWidth', 2)
    ylabel('Phase [deg]', 'FontSize', 14);
    xlabel('Freq [Hz]', 'FontSize', 14);
    axis([min(f) max(f) 90*floor(min(angle(G)*180/pi)/90) 90*ceil(max(angle(G)*180/pi)/90)])
%    set(gca,'YTick',[-1800:90:1800])
    grid on
    if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
    end

    
%%  % FSS Close Loop Response Response
    hdl= figure(103)
    a = SYS(21,19);
    b = SYS(22,19);
    G = mybodesys(b/a,f);
    %
    subplot(211)
    semilogx(f,20*log10(abs(G)), 'LineWidth', 2)
    tt= title('FSS Open Loop TF (in 19, out 22/out 21)', 'FontSize', 16);
    ylabel('Mag [dB]', 'FontSize', 14);
    axis([min(f) max(f) floor(min(log10(abs(G))))*20 ceil(max(log10(abs(G))))*20])
    grid on
    subplot(212)
    semilogx(f,angle(G)*180/pi, 'LineWidth', 2)
    ylabel('Phase [deg]', 'FontSize', 14);
    xlabel('Freq [Hz]', 'FontSize', 14);
    axis([min(f) max(f) 90*floor(min(angle(G)*180/pi)/90) 90*ceil(max(angle(G)*180/pi)/90)])
%    set(gca,'YTick',[-1800:90:1800])
    grid on
    if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
    end

%%  % FSS Supression Response
    hdl= figure(104)
    G = mybodesys(SYS(21,19),f);
    %
    subplot(211)
    semilogx(f,20*log10(abs(G)), 'LineWidth', 2)
    tt= title('FSS Suppression Response (in 19, out 21)', 'FontSize', 16);
    ylabel('Mag [dB]', 'FontSize', 14);
    axis([min(f) max(f) floor(min(log10(abs(G))))*20 ceil(max(log10(abs(G))))*20])
    grid on
    subplot(212)
    semilogx(f,angle(G)*180/pi, 'LineWidth', 2)
    ylabel('Phase [deg]', 'FontSize', 14);
    xlabel('Freq [Hz]', 'FontSize', 14);
    axis([min(f) max(f) 90*floor(min(angle(G)*180/pi)/90) 90*ceil(max(angle(G)*180/pi)/90)])
%    set(gca,'YTick',[-1800:90:1800])
    grid on
    if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
    end
    
end         % end FSS



if TTFSS
%% TTFSS Feedback of the laser to the Heterodyen Signal
    %
    % Local laser to TTFSS error signal
    hdl= figure(1)
    a = SYS(5,8);
    b = SYS(3,8);
    G = mybodesys(b/a,f);
    %
    subplot(211)
    semilogx(f,20*log10(abs(G)), 'LineWidth', 2)
    tt= title('Local laser -to- TTFSS error signal TF (in 8, out 3/out 5)', 'FontSize', 16);
    ylabel('Mag [dB]', 'FontSize', 14);
    axis([min(f) max(f) floor(min(log10(abs(G))))*20 ceil(max(log10(abs(G))))*20])
    grid on
    subplot(212)
    semilogx(f,angle(G)*180/pi, 'LineWidth', 2)
    ylabel('Phase [deg]', 'FontSize', 14);
    xlabel('Freq [Hz]', 'FontSize', 14);
    axis([min(f) max(f) 90*floor(min(angle(G)*180/pi)/90) 90*ceil(max(angle(G)*180/pi)/90)])
%    set(gca,'YTick',[-1800:90:1800])
    grid on
    if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
    end


    % TTFSS Open Loop response
    hdl= figure(2)
    a = SYS(6,7);
    b = SYS(5,7);
    G = mybodesys(b/a,f);
    %
    subplot(211)
    semilogx(f,20*log10(abs(G)), 'LineWidth', 2)
    tt= title('TTFSS Controller TF (in 7, out 5/out 6)','FontSize', 16);
    ylabel('Mag [dB]', 'FontSize', 14);
    axis([min(f) max(f) floor(min(log10(abs(G))))*20 ceil(max(log10(abs(G))))*20])
    grid on
    subplot(212)
    semilogx(f,angle(G)*180/pi, 'LineWidth', 2)
    ylabel('Phase [deg]', 'FontSize', 14);
    xlabel('Freq [Hz]', 'FontSize', 14);
    axis([min(f) max(f) 90*floor(min(angle(G)*180/pi)/90) 90*ceil(max(angle(G)*180/pi)/90)])
%    set(gca,'YTick',[-1800:90:1800])
    grid on
    if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
    end

    
    % TTFSS Close Loop Response Response
    hdl= figure(3)
    a = SYS(6,7);
    b = SYS(3,7);
    G = mybodesys(b/a,f);
    %
    subplot(211)
    semilogx(f,20*log10(abs(G)), 'LineWidth', 2)
    tt= title('TTFSS Open Loop TF (in 7, out 3/ out 6)', 'FontSize', 16);
    ylabel('Mag [dB]', 'FontSize', 14);
    axis([min(f) max(f) floor(min(log10(abs(G))))*20 ceil(max(log10(abs(G))))*20])
    grid on
    subplot(212)
    semilogx(f,angle(G)*180/pi, 'LineWidth', 2)
    ylabel('Phase [deg]', 'FontSize', 14);
    xlabel('Freq [Hz]', 'FontSize', 14);
    axis([min(f) max(f) 90*floor(min(angle(G)*180/pi)/90) 90*ceil(max(angle(G)*180/pi)/90)])
%    set(gca,'YTick',[-1800:90:1800])
    grid on
    if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
    end

    % TTFSS Supression Response
    hdl= figure(4)
    G = mybodesys(SYS(6,7),f);
    %
    subplot(211)
    semilogx(f,20*log10(abs(G)), 'LineWidth', 2)
    tt= title('TTFSS Suppression Response (in 7, out 6)', 'FontSize', 16);
    ylabel('Mag [dB]', 'FontSize', 14);
    axis([min(f) max(f) floor(min(log10(abs(G))))*20 ceil(max(log10(abs(G))))*20])
    grid on
    subplot(212)
    semilogx(f,angle(G)*180/pi, 'LineWidth', 2)
    ylabel('Phase [deg]', 'FontSize', 14);
    xlabel('Freq [Hz]', 'FontSize', 14);
    axis([min(f) max(f) 90*floor(min(angle(G)*180/pi)/90) 90*ceil(max(angle(G)*180/pi)/90)])
%    set(gca,'YTick',[-1800:90:1800])
    grid on
    if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
    end
    
end         % end TTFSS

if PDH
%% PDH Feedback of the laser frequency to the arm cavity, via the VCO,EX
    %    
    % Local VCO,EX to PDH error signal
    hdl = figure(5)
    a = SYS(8,9);
    b = SYS(10,9);
    G = mybodesys(b/a,f);
    %
    subplot(211)
    semilogx(f,20*log10(abs(G)), 'LineWidth', 2)
    tt= title('VCO,EX -to- PDH Error Signal TF (in 9, out 10/out 8)', 'FontSize', 16);
    ylabel('Mag [dB]', 'FontSize', 14);
    axis([min(f) max(f) floor(min(log10(abs(G))))*20 ceil(max(log10(abs(G))))*20])
    grid on
    subplot(212)
    semilogx(f,angle(G)*180/pi, 'LineWidth', 2)
    ylabel('Phase [deg]', 'FontSize', 14);
    xlabel('Freq [Hz]', 'FontSize', 14);
    axis([min(f) max(f) 90*floor(min(angle(G)*180/pi)/90) 90*ceil(max(angle(G)*180/pi)/90)])
%    set(gca,'YTick',[-1800:90:1800])
    grid on
    if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
    end
    
    
    % PDH Controller response
    hdl= figure(6)
    a = SYS(9,10);
    b = SYS(8,10);
    G = mybodesys(b/a,f);
    %
    subplot(211)
    semilogx(f,20*log10(abs(G)), 'LineWidth', 2)
    tt= title('PDH Controller TF (in 10, out 8/out 9)', 'FontSize', 16);
    ylabel('Mag [dB]', 'FontSize', 14);
    axis([min(f) max(f) floor(min(log10(abs(G))))*20 ceil(max(log10(abs(G))))*20])
    grid on
    subplot(212)
    semilogx(f,angle(G)*180/pi, 'LineWidth', 2)
    ylabel('Phase [deg]', 'FontSize', 14);
    xlabel('Freq [Hz]', 'FontSize', 14);
    axis([min(f) max(f) 90*floor(min(angle(G)*180/pi)/90) 90*ceil(max(angle(G)*180/pi)/90)])
%    set(gca,'YTick',[-1800:90:1800])
    grid on
    if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
    end

    
    % PDH Close Loop Response Response
    hdl = figure(7)
    a = SYS(9,10);
    b = SYS(10,10);
    G = mybodesys(b/a,f);
    %
    subplot(211)
    semilogx(f,20*log10(abs(G)), 'LineWidth', 2)
    tt= title('PDH Open Loop TF (in 10, out 10/out 9)', 'FontSize', 16);
    ylabel('Mag [dB]', 'FontSize', 14);
    axis([min(f) max(f) floor(min(log10(abs(G))))*20 ceil(max(log10(abs(G))))*20])
    grid on
    subplot(212)
    semilogx(f,angle(G)*180/pi, 'LineWidth', 2)
    ylabel('Phase [deg]', 'FontSize', 14);
    xlabel('Freq [Hz]', 'FontSize', 14);
    axis([min(f) max(f) 90*floor(min(angle(G)*180/pi)/90) 90*ceil(max(angle(G)*180/pi)/90)])
%    set(gca,'YTick',[-1800:90:1800])
    grid on
    if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
    end

    % PDH Supression Response
    hdl = figure(8)
    G = mybodesys(SYS(9,10),f);
    %
    subplot(211)
    semilogx(f,20*log10(abs(G)), 'LineWidth', 2)
    tt= title('PDH Suppression Response (in 10, out 9)', 'FontSize', 16);
    ylabel('Mag [dB]', 'FontSize', 14);
    axis([min(f) max(f) floor(min(log10(abs(G))))*20 ceil(max(log10(abs(G))))*20])
    grid on
    subplot(212)
    semilogx(f,angle(G)*180/pi, 'LineWidth', 2)
    ylabel('Phase [deg]', 'FontSize', 14);
    xlabel('Freq [Hz]', 'FontSize', 14);
    axis([min(f) max(f) 90*floor(min(angle(G)*180/pi)/90) 90*ceil(max(angle(G)*180/pi)/90)])
%    set(gca,'YTick',[-1800:90:1800])
    grid on
    if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
    end

end         % end PDH

if COMM
%% Common Mode Feedback of the PSL frequency to the common mode arm cavity
%% length fluctuations, via the VCO,C
    %
    % Vertex VCO,PSL to Common Mode error signal
    hdl = figure(9)
    a = SYS(12,12);
    b = SYS(13,12);
    G = mybodesys(b/a,f);
    %
    subplot(211)
    semilogx(f,20*log10(abs(G)), 'LineWidth', 2)
    tt= title('VCO,PSL -to- Common Mode Error Signal TF (in 12, out 13/out 13)', 'FontSize', 16);
    ylabel('Mag [dB]', 'FontSize', 14);
    axis([min(f) max(f) floor(min(log10(abs(G))))*20 ceil(max(log10(abs(G))))*20])
    grid on
    subplot(212)
    semilogx(f,angle(G)*180/pi, 'LineWidth', 2)
    ylabel('Phase [deg]', 'FontSize', 14);
    xlabel('Freq [Hz]', 'FontSize', 14);
    axis([min(f) max(f) 90*floor(min(angle(G)*180/pi)/90) 90*ceil(max(angle(G)*180/pi)/90)])
%    set(gca,'YTick',[-1800:90:1800])
    grid on
    if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
    end

    % Common Mode Servo Open Loop response
    hdl= figure(10)
    a = SYS(11,11);
    b = SYS(12,11);
    G = mybodesys(b/a,f);
    %
    subplot(211)
    semilogx(f,20*log10(abs(G)), 'LineWidth', 2)
    tt= title('Common Mode Controller TF (in 11, out 12/ out 11)', 'FontSize',16);
    ylabel('Mag [dB]', 'FontSize', 14)
    axis([min(f) max(f) floor(min(log10(abs(G))))*20 ceil(max(log10(abs(G))))*20])
    grid on
    subplot(212)
    semilogx(f,angle(G)*180/pi, 'LineWidth', 2)
    ylabel('Phase [deg]', 'FontSize', 14);
    xlabel('Freq [Hz]', 'FontSize', 14);
    axis([min(f) max(f) 90*floor(min(angle(G)*180/pi)/90) 90*ceil(max(angle(G)*180/pi)/90)])
%    set(gca,'YTick',[-1800:90:1800])
    grid on
    if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
    end
    
    
    % Common Mode Close Loop Response
    hdl = figure(11)
    a = SYS(11,11);
    b = SYS(13,11);
    G = mybodesys(b/a,f);
    %
    subplot(211)
    semilogx(f,20*log10(abs(G)), 'LineWidth', 2)
    tt= title('Common Mode Open Loop TF (in 11, out 13/out 11)', 'FontSize', 16);
    ylabel('Mag [dB]', 'FontSize', 14);
    axis([min(f) max(f) floor(min(log10(abs(G))))*20 ceil(max(log10(abs(G))))*20])
    grid on
    subplot(212)
    semilogx(f,angle(G)*180/pi, 'LineWidth', 2)
    ylabel('Phase [deg]', 'FontSize', 14);
    xlabel('Freq [Hz]', 'FontSize', 14);
    axis([min(f) max(f) 90*floor(min(angle(G)*180/pi)/90) 90*ceil(max(angle(G)*180/pi)/90)])
%    set(gca,'YTick',[-1800:90:1800])
    grid on
    if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
    end
%%
    % Common Mode Supression Response
    hdl = figure(12)
    G = mybodesys(SYS(11,11),f);
    %
    subplot(211)
    semilogx(f,20*log10(abs(G)), 'LineWidth', 2)
    tt= title('Common Mode Suppression Response (in 11, out 11)', 'FontSize', 16);
    ylabel('Mag [dB]', 'FontSize', 14);
    axis([min(f) max(f) floor(min(log10(abs(G))))*20 ceil(max(log10(abs(G))))*20])
    grid on
    subplot(212)
    semilogx(f,angle(G)*180/pi, 'LineWidth', 2)
    ylabel('Phase [deg]', 'FontSize', 14);
    xlabel('Freq [Hz]', 'FontSize', 14);
    axis([min(f) max(f) 90*floor(min(angle(G)*180/pi)/90) 90*ceil(max(angle(G)*180/pi)/90)])
%    set(gca,'YTick',[-1800:90:1800])
    grid on
    if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
    end
end         % end COMM

if DIFF
%% Differential Mode Feedback to both the ETMs (out of phase). This
%% requires the Quad response and its servo, for now I have a single
%% pendulum replacing the Quad...
    %    
    % Diff Mode Servo input to Differential Mode error signal
    hdl = figure(13)
    a = SYS(23,20);
    b = SYS(14,20);
    G = mybodesys(b/a,f);
    %
    subplot(211)
    semilogx(f,20*log10(abs(G)), 'LineWidth', 2)
    tt= title('Diff Servo Ouput -to- Diff Mode Error Signal TF (in 20, out 14/out 23)', 'FontSize', 16);
    ylabel('Mag [dB]', 'FontSize', 14);
%    axis([min(f) max(f) floor(min(log10(abs(G))))*20 ceil(max(log10(abs(G))))*20])
    axis([min(f) max(f) -50 250])
    grid on
    subplot(212)
    semilogx(f,angle(G)*180/pi, 'LineWidth', 2)
    ylabel('Phase [deg]', 'FontSize', 14);
    xlabel('Freq [Hz]', 'FontSize', 14);
    axis([min(f) max(f) 90*floor(min(angle(G)*180/pi)/90) 90*ceil(max(angle(G)*180/pi)/90)])
%    set(gca,'YTick',[-1800:90:1800])
    grid on
    if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
    end

    
    % Diff Mode Controller response
    hdl= figure(14)
    a = SYS(17,15);
    b = SYS(23,15);
    G = mybodesys(b/a,f);
    %
    subplot(211)
    semilogx(f,20*log10(abs(G)), 'LineWidth', 2)
    tt= title('Differential Mode Controller TF (in 15, out 23/out 17)', 'FontSize',16);
    ylabel('Mag [dB]', 'FontSize', 14)
    axis([min(f) max(f) floor(min(log10(abs(G))))*20 ceil(max(log10(abs(G))))*20])
    grid on
    subplot(212)
    semilogx(f,angle(G)*180/pi, 'LineWidth', 2)
    ylabel('Phase [deg]', 'FontSize', 14);
    xlabel('Freq [Hz]', 'FontSize', 14);
    axis([min(f) max(f) 90*floor(min(angle(G)*180/pi)/90) 90*ceil(max(angle(G)*180/pi)/90)])
%    set(gca,'YTick',[-1800:90:1800])
    grid on
    if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
    end
    
    % Differentail Mode Close Loop Response
    hdl = figure(15)
    a = SYS(17,15);
    b = SYS(18,15);
    G = mybodesys(b/a,f);
    %
    subplot(211)
    semilogx(f,20*log10(abs(G)), f, 20*log10(50./f), 'LineWidth', 2)
    tt= title('Differential Mode Open Loop TF (in 15, out 18/out 17)', 'FontSize', 16);
    ylabel('Mag [dB]', 'FontSize', 14);
%    axis([min(f) max(f) floor(min(log10(abs(G))))*20 ceil(max(log10(abs(G))))*20])
    axis([min(f) max(f) -100 100])
    grid on
    subplot(212)
    semilogx(f,angle(G)*180/pi, 'LineWidth', 2)
    ylabel('Phase [deg]', 'FontSize', 14);
    xlabel('Freq [Hz]', 'FontSize', 14);
    axis([min(f) max(f) 90*floor(min(angle(G)*180/pi)/90) 90*ceil(max(angle(G)*180/pi)/90)])
%    set(gca,'YTick',[-1800:90:1800])
    grid on
    if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
    end

    % Differential Mode Supression Response
    hdl = figure(16)
    G = mybodesys(SYS(17,15),f);
    %
    subplot(211)
    semilogx(f,20*log10(abs(G)), 'LineWidth', 2)
    tt= title('Differential Mode Suppression Response (in 15, out 17)', 'FontSize', 16);
    ylabel('Mag [dB]', 'FontSize', 14);
    axis([min(f) max(f) floor(min(log10(abs(G))))*20 ceil(max(log10(abs(G))))*20])
    grid on
    subplot(212)
    semilogx(f,angle(G)*180/pi, 'LineWidth', 2)
    ylabel('Phase [deg]', 'FontSize', 14);
    xlabel('Freq [Hz]', 'FontSize', 14);
    axis([min(f) max(f) 90*floor(min(angle(G)*180/pi)/90) 90*ceil(max(angle(G)*180/pi)/90)])
%    set(gca,'YTick',[-1800:90:1800])
    grid on
    if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
    end
    
end         % end DIFF
%%
run = 1;
if run
    a = mybodesys(SYS(36,6), f); % Ref Cav noise -to- PDH Freq Shift
    b = mybodesys(SYS(4,6), f); % Ref Cav noise -to- Local Laser Freq Shift (TTFSS)
    c = mybodesys(SYS(34,6), f); % Ref Cav noise -to- COMM Freq Shift
    d = mybodesys(SYS(35,6), f); % Ref Cav noise -to- Diff Freq Shift
    e = mybodesys(SYS(25,6), f); % Ref Cav noise -to- X-arm freq Shift
    
    PDHfreq = a .* NoisyRefCav;
    TTFSSfreq = b  .* NoisyRefCav;
    COMMfreq = c .* NoisyRefCav;
    DIFFfreq = d .* NoisyRefCav;
    XarmFreq = e .* NoisyRefCav;

    
    G = [PDHfreq TTFSSfreq COMMfreq DIFFfreq];

    hdl = figure(20)
    loglog(f,(abs(G)), ...
        f, abs(NoisyRefCav), ...
        f, rms(f, abs(NoisyRefCav)'), '--','LineWidth', 2)
    legend('PDH freq, out(36)', 'TTFSS freq, out(4)', 'COMM freq, out(34)', 'DIFF freq, out(35)', ...
        'Noisy Reference Cavity, in(6)', 'Noisy Reference Cavity, rms');
    tt= title('ALS System response with noisy Ref Cavity.', 'FontSize', 18);
    ylabel('Frequency Noise [Hz/rtHz]', 'FontSize', 16);
    axis([min(f) max(f) 1e-7 1e3]);
    set(gca, 'FontSize', 14);
    grid on
    xlabel('Freq [Hz]', 'FontSize', 16);
    grid on
%    if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
%    end
    
end

%%
    hdl = figure(21)
    loglog(f,(abs(G).*L_arm/nu_GRN), ...
        f, abs(NoisyRefCav).*L_arm/nu_GRN, ...
        f, rms(f, abs(NoisyRefCav)').*L_arm/nu_GRN, '--','LineWidth', 2)
    legend('PDH freq, out(36)', 'TTFSS freq, out(4)', 'COMM freq, out(34)', 'DIFF freq, out(35)', ...
        'Noisy Reference Cavity, in(6)', 'Noisy Reference Cavity, rms');
    tt= title('ALS System response with noisy Ref Cavity.', 'FontSize', 18);
    ylabel('Displacement Noise [m/rtHz]', 'FontSize', 16);
    axis([min(f) max(f) 1e-20 1e-8]);
    set(gca, 'FontSize', 14);
    grid on
    xlabel('Freq [Hz]', 'FontSize', 16);
    grid on
    
%% Addition for the ALS publication

SEI_F = [0.01 0.03 0.1 0.2 0.5 1 10 30 100];
SEI_X = [3e-6 1e-6 2e-7 2e-7 8e-10 1e-11 3e-13 3e-14 5e-15];
sei_model_0612 = 10.^(interp1(SEI_F,log10(SEI_X),f,'cubic',-14));

a = mybodesys(SYS(24,30), f);
%
residual = abs(a).* sei_model_0612';
residual_rms = rms(f, residual');
%%
% figure(99); loglog(f, residual);
% OpenloopResidual = residual;
% OpenloopResidual_rms = residual_rms;
% save('openloopresiduals.mat', 'OpenloopResidual', 'OpenloopResidual_rms');

%%
hdl = figure(22)
close(hdl);
load('openloopresiduals');
res = [residual tmOpenloopRes]';
res_rms = [residual_rms' tmOpenloopRes_rms']';

vel = abs(tmOpenloopRes' .* i*2*pi*f);
vel_rms = rms(f, vel);

hdl = figure(22)
%
[AX, H1, H2] = plotyy(f, res, ...
                      f, res_rms, 'loglog');
set(get(AX(1), 'YLabel'), 'String', ['Displacement Noise Spectrum [m/rtHz]'], 'FontSize', 20);
set(get(AX(2), 'YLabel'), 'String', ['Displacement Noise RMS [m]'], 'FontSize', 20);
xlabel('Freq [Hz]', 'FontSize', 20);
grid on
legend('Closeloop Test mass','Openloop Test mass','Closeloop Test mass RMS', 'Openloop Test mass RMS');
set(H1, 'LineStyle', '-', 'LineWidth', 2);
set(H2, 'LineStyle', '--', 'LineWidth', 2);
set(AX(1), 'ylim', [1e-15 1e-5]);
set(AX(1), 'ytickmode', 'auto');
linkaxes([AX(1), AX(2)], 'y');
set(AX(1), 'ytick', [1e-15 1e-14 1e-13 1e-12 1e-11 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5]);
set(AX(2), 'ytick', [1e-14 1e-13 1e-12 1e-11 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5]);
set(AX(1), 'FontSize', 18);
set(AX(2), 'FontSize', 18);
%
% hdl = loglog(f, res, f, res_rms, ...)
%              'LineWidth', 2)
% legend('PDH freq, out(36)', 'TTFSS freq, out(4)', 'COMM freq, out(34)', 'DIFF freq, out(35)', ...
%        'Noisy Reference Cavity, in(6)', 'Noisy Reference Cavity, rms');
%    tt= title('ALS System response with noisy Ref Cavity.', 'FontSize', 18);
%    ylabel('Displacement Noise [m/rtHz]', 'FontSize', 16);
%    axis([min(f) max(f) 1e-5 1e0]);
%    set(gca, 'FontSize', 14);
%    grid on
%    xlabel('Freq [Hz]', 'FontSize', 16);
%    grid on
%     if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
        print('-dpdf', [save_figure_dir , 'LockingStrategy05_TM_', num2str(hdl,'%.3d'), '.pdf']);
%     end


%% X-arm Beam line motion to Common Mode Frequency Fluctuations
% input: 30
% output: 34
% sys(output, input);

a = mybodesys(SYS(34,30), f);
comm_residual = abs(a).* L_arm .* sei_model_0612' ./ nu_GRN;
comm_residual_rms = rms(f, comm_residual');

% X-arm Beam line motion to Differential Mode Frequency Fluctuations
% input: 30
% output: 35
% sys(output, input);

a = mybodesys(SYS(35,30), f);
diff_residual = abs(a).* L_arm .* sei_model_0612' ./ nu_GRN;
diff_residual_rms = rms(f, diff_residual');

res = [comm_residual diff_residual]';
res_rms = [comm_residual_rms' diff_residual_rms']';

hdl = figure(23)
[AX, H1, H2] = plotyy(f, res, ...
                      f, res_rms, 'loglog');
set(get(AX(1), 'YLabel'), 'String', ['Displacement Noise Spectrum [m/rtHz]'], 'FontSize', 20);
set(get(AX(2), 'YLabel'), 'String', ['Displacement Noise RMS [m]'], 'FontSize', 20);
xlabel('Freq [Hz]', 'FontSize', 20);
grid on
legend('Common Mode','Differential Mode','Common Mode RMS', 'Differential Mode RMS');
set(H1, 'LineStyle', '-', 'LineWidth', 2);
set(H2, 'LineStyle', '--', 'LineWidth', 2);
set(AX(1), 'ylim', [1e-21 1e-13]);
set(AX(1), 'ytickmode', 'auto');
linkaxes([AX(1), AX(2)], 'y');
set(AX(1), 'ytick', [1e-21 1e-20 1e-19 1e-18 1e-17 1e-16 1e-15 1e-14 1e-13]);
set(AX(2), 'ytick', [1e-20 1e-19 1e-18 1e-17 1e-16 1e-15 1e-14 1e-13]);
set(AX(1), 'FontSize', 18);
set(AX(2), 'FontSize', 18);

%     if save_figure == 1
        orient(hdl, 'landscape');
        print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
        print('-dpdf', [save_figure_dir , 'LockingStrategy05_CommDiff_', num2str(hdl,'%.3d'), '.pdf']);
%     end

%% PSL Doubling phase fluctuations due to crystal temperature stabilities
% EXC input: 31 (at the frequency location, not phase)
% IN1 out: 37
% IN2 out: 38
%
% sys(output, input)
%

a = mybodesys(SYS(36,31), f); % PSL-SHG,PHS -to- PDH Freq Shift
b = mybodesys(SYS(4,31), f); % PSL-SHG,PHS -to- Local Laser Freq Shift (TTFSS)
c = mybodesys(SYS(34,31), f); % PSL-SHG,PHS -to- COMM Freq Shift
d = mybodesys(SYS(35,31), f); % PSL-SHG,PHS -to- Diff Freq Shift
e = mybodesys(SYS(25,31), f); % PSL-SHG,PHS -to- X-arm freq Shift


SHGPhaseNoise = ones(size(f)) .* 2*pi;

SHGFreqNoise = SHGPhaseNoise ./ f;
SHG_Equiv_DisplNoise = SHGFreqNoise' .* L_arm / nu_IR;

PDHfreq = a .* SHG_Equiv_DisplNoise;
TTFSSfreq = b  .* SHG_Equiv_DisplNoise;
COMMfreq = c .* SHG_Equiv_DisplNoise;
DIFFfreq = d .* SHG_Equiv_DisplNoise;
XarmFreq = e .* SHG_Equiv_DisplNoise;


G = [PDHfreq TTFSSfreq COMMfreq DIFFfreq];

hdl = figure(24)
loglog(f,(abs(G)), ...
    f, abs(SHG_Equiv_DisplNoise), ...
    f, rms(f, abs(SHG_Equiv_DisplNoise)'), '--','LineWidth', 2)
legend('PDH freq, out(36)', 'TTFSS freq, out(4)', 'COMM freq, out(34)', 'DIFF freq, out(35)', ...
    'PSL SHG Freq Noise, in(31)', 'PSL SHG Freq Noise, rms', ...
    'Location', 'NorthWest');
tt= title({'ALS System response from PSL SHG Frequency Drift' '(due to crystal temperature fluctuations, ~10 mrad)'}, 'FontSize', 18);
ylabel('Displacement Noise [m/rtHz]', 'FontSize', 16);
axis([min(f) max(f) 1e-20 1e-14]);
set(gca, 'FontSize', 14);
grid on
xlabel('Freq [Hz]', 'FontSize', 16);
grid on
%    if save_figure == 1
    orient(hdl, 'landscape');
    print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
%    end
