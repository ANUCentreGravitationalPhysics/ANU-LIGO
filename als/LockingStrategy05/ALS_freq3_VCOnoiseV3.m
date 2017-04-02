% It is assumed that ALS_lockingv3 has been run.
%

save_figure = 0;
% save_figure_dir = 'sim/';

f = logspace(-1,4.4, 1e3);

Pt = 10e-3; % power of carrier * power of sidebands
PDHslope = 16*sqrt(Pt)*ArmCavFIN/lambda_GRN;

%% Load VCO Noise Model
vco = load('T0900451_Model_FirstSecond_Stage.txt');
f_vco = vco(:,1);
firststage_dBc = vco(:,2);
secondstage_dBc = vco(:,3);

firststage_Hz = sqrt(2 .* 10.^(firststage_dBc/10) ) .* f_vco;
secondstage_Hz = sqrt(2 .* 10.^(secondstage_dBc/10) ) .* f_vco;

nu_GRN = cc / lambda_GRN;
firststage_m = firststage_Hz * L_arm / nu_GRN;
secondstage_m = secondstage_Hz * L_arm / nu_GRN;

hdl = figure(999)
loglog(f_vco, firststage_m, f_vco, secondstage_m, 'LineWidth', 2);
grid on;
set(gca, 'FontSize', 14);
legend('Model (First Stage)', 'Model (Second Stage)');
title({'VCO Noise (SSB)' 'See T0900451'}, 'FontSize', 16);
xlabel('Frequency [Hz]', 'FontSize', 16);
ylabel('Phase noise /L_{arm} [m/rtHz]', 'FontSize', 16);
if save_figure == 1
    orient(hdl, 'landscape');
%    print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), ' ',get(tt,'string'), '.pdf']);
    print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
end


%% VCO input to ETMx equivalent frequency shift 
hdl = figure(1000)
a = mybodesys(SYS(25,26), f);  % PSL,VCO to ETMx equivalent Freq Shift, Hz/Hz
b = mybodesys(SYS(25,27), f);  % f_AOM to ETMx, Hz/Hz
c = mybodesys(SYS(25,25), f);  % VCO,EX to ETMx
d = mybodesys(SYS(25,28), f); % COMM,VCO to ETMx
e = mybodesys(SYS(25,29), f); % DIFF,VCO to ETMx
%
subplot(211);
set(gca, 'FontSize', 14);
semilogx(f,20*log10(abs(a)), ...
    f,20*log10(abs(b)), ...
    f,20*log10(abs(c)), ...
    f,20*log10(abs(d)), ...
    f,20*log10(abs(e)), 'LineWidth', 2)
tt = title('VCOs -to- ETMx Equivalent Freq. Shift (out 25)','FontSize', 16);
legend('PSL,VCO (in 26)', 'AOM,VCO (in 27)', 'EX,VCO (in 25)', 'COMM,VCO (in 28)', 'DIFF,VCO (in 29)')
ylabel('Mag [dB]', 'FontSize', 14);
axis([min(f) max(f) -240 60]);
set(gca, 'FontSize', 14);
grid on
subplot(212);
semilogx(f, angle(a)*180/pi, ...
         f, angle(b)*180/pi, ...
         f, angle(c)*180/pi, ...
         f, angle(d)*180/pi, ...
         f, angle(e)*180/pi, 'LineWidth', 2)
ylabel('Phase [deg]', 'FontSize', 14);
xlabel('Freq [Hz]', 'FontSize', 14);
axis([min(f) max(f) 90*floor(min(angle(b)*180/pi)/90) 90*ceil(max(angle(b)*180/pi)/90)])
set(gca,'YTick',[-1800:90:1800])
set(gca, 'FontSize', 14);

grid on
if save_figure == 1
    orient(hdl, 'landscape');
%    print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), ' ',get(tt,'string'), '.pdf']);
    print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
end

% VCO input to Common Mode Error Signal
hdl = figure(1001);
a = mybodesys(SYS(34,26), f);  % PSL,VCO (Volts) to Comm Freq Shift
b = mybodesys(SYS(34,27), f);  % f_AOM to Comm Freq Shift
c = mybodesys(SYS(34,25), f);  % VCO,EX to Comm Freq Shift
d = mybodesys(SYS(34,28), f); % COMM,VCO to Comm Freq Shift
e = mybodesys(SYS(34,29), f); % DIFF,VCO to Comm Freq Shift
%
subplot(211);
semilogx(f,20*log10(abs(a)), ...
    f,20*log10(abs(b)), ...
    f,20*log10(abs(c)), ...
    f,20*log10(abs(d)), ...
    f,20*log10(abs(e)), 'LineWidth', 2)
tt = title({'VCOs -to- Common Mode Freq. Shift (out 34)', 'Freq shift between PSL and X-arm Laser'},'FontSize', 16);
legend('PSL,VCO (in 26)', 'AOM,VCO (in 27)', 'EX,VCO (in 25)', 'COMM,VCO (in 28)', 'DIFF,VCO (in 29)', ...
        'Location', 'SouthEast')
ylabel('Mag [dB]', 'FontSize', 14);
axis([min(f) max(f) -220 20])
grid on
subplot(212);
semilogx(f, angle(a)*180/pi, ...
         f, angle(b)*180/pi, ...
         f, angle(c)*180/pi, ...
         f, angle(d)*180/pi, ...
         f, angle(e)*180/pi, 'LineWidth', 2)
ylabel('Phase [deg]', 'FontSize', 14);
xlabel('Freq [Hz]', 'FontSize', 14);
axis([min(f) max(f) 90*floor(min(angle(b)*180/pi)/90) 90*ceil(max(angle(b)*180/pi)/90)])
    set(gca,'YTick',[-1800:90:1800])
grid on
if save_figure == 1
    orient(hdl, 'landscape');
%    print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), ' ',get(tt,'string'), '.pdf']);
    print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
end

% VCO input to Differential Mode Error Signal
hdl = figure(1002);
a = mybodesys(SYS(35,26), f);  % PSL,VCO (Volts) to Diff Freq Shift
b = mybodesys(SYS(35,27), f);  % f_AOM to Diff Freq Shift
c = mybodesys(SYS(35,25), f);  % VCO,EX to Diff Freq Shift
d = mybodesys(SYS(35,28), f); % COMM,VCO to Diff Freq Shift
e = mybodesys(SYS(35,29), f); % DIFF,VCO to Diff Freq Shift
%
subplot(211);
semilogx(f,20*log10(abs(a)), ...
    f,20*log10(abs(b)), ...
    f,20*log10(abs(c)), ...
    f,20*log10(abs(d)), ...
    f,20*log10(abs(e)), 'LineWidth', 2)
tt = title({'VCOs -to- Differential Mode Freq. Shift (out 35)', 'Freq shift between X-arm and Y-arm Lasers'},'FontSize', 16);
legend('PSL,VCO (in 26)', 'AOM,VCO (in 27)', 'EX,VCO (in 25)', 'COMM,VCO (in 28)', 'DIFF,VCO (in 29)', ...
        'Location', 'SouthEast')
ylabel('Mag [dB]', 'FontSize', 14);
axis([min(f) max(f) -350 20])
grid on
subplot(212);
semilogx(f, angle(a)*180/pi, ...
         f, angle(b)*180/pi, ...
         f, angle(c)*180/pi, ...
         f, angle(d)*180/pi, ...
         f, angle(e)*180/pi, 'LineWidth', 2)
ylabel('Phase [deg]', 'FontSize', 14);
xlabel('Freq [Hz]', 'FontSize', 14);
axis([min(f) max(f) 90*floor(min(angle(b)*180/pi)/90) 90*ceil(max(angle(b)*180/pi)/90)])
    set(gca,'YTick',[-1800:90:1800])
grid on
if save_figure == 1
    orient(hdl, 'landscape');
%    print('-dpdf', [save_figure_dir,num2str(hdl,'%.3d'), ' ',get(tt,'string'), '.pdf']);
    print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
end

%% VCO input to Local Laser Freq
hdl = figure(1003);
a = mybodesys(SYS(4,26), f);  % PSL,VCO (Volts) to TTFSS Local Laser Freq Shift
b = mybodesys(SYS(4,27), f);  % f_AOM to TTFSS Local Laser Freq Shift
c = mybodesys(SYS(4,25), f);  % VCO,EX to TTFSS Local Laser Freq Shift
d = mybodesys(SYS(4,28), f); % COMM,VCO to TTFSS Local Laser Freq Shift
e = mybodesys(SYS(4,29), f); % DIFF,VCO to TTFSS Local Laser Freq Shift
%
subplot(211);
semilogx(f,20*log10(abs(a)), ...
    f,20*log10(abs(b)), ...
    f,20*log10(abs(c)), ...
    f,20*log10(abs(d)), ...
    f,20*log10(abs(e)), 'LineWidth', 2)
tt = title({'VCOs -to- Local Laser Freq (out 4)'},'FontSize', 16);
legend('PSL,VCO (in 26)', 'AOM,VCO (in 27)', 'EX,VCO (in 25)', 'COMM,VCO (in 28)', 'DIFF,VCO (in 29)', ...
        'Location', 'SouthEast')
ylabel('Mag [dB]', 'FontSize', 14);
axis([min(f) max(f) -300 40])
grid on
subplot(212);
semilogx(f, angle(a)*180/pi, ...
         f, angle(b)*180/pi, ...
         f, angle(c)*180/pi, ...
         f, angle(d)*180/pi, ...
         f, angle(e)*180/pi, 'LineWidth', 2)
ylabel('Phase [deg]', 'FontSize', 14);
xlabel('Freq [Hz]', 'FontSize', 14);
axis([min(f) max(f) 90*floor(min(angle(b)*180/pi)/90) 90*ceil(max(angle(b)*180/pi)/90)])
    set(gca,'YTick',[-1800:90:1800])
grid on
if save_figure == 1
    orient(hdl, 'landscape');
%    print('-dpdf', [save_figure_dir,num2str(hdl,'%.3d'), ' ',get(tt,'string'), '.pdf']);
    print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
end

%% VCO input to PDH Error Signal
hdl = figure(1004);
a = mybodesys(SYS(36,26), f);  % PSL,VCO (Volts) to PDH Freq Shift
b = mybodesys(SYS(36,27), f);  % f_AOM to PDH Freq Shift
c = mybodesys(SYS(36,25), f);  % VCO,EX to PDH Freq Shift
d = mybodesys(SYS(36,28), f); % COMM,VCO to PDH Freq Shift
e = mybodesys(SYS(36,29), f); % DIFF,VCO to PDH Freq Shift
%
subplot(211);
semilogx(f,20*log10(abs(a)), ...
    f,20*log10(abs(b)), ...
    f,20*log10(abs(c)), ...
    f,20*log10(abs(d)), ...
    f,20*log10(abs(e)), 'LineWidth', 2)
tt = title({'VCOs -to- PDH Freq Shift (out 10)', 'Freq shift away from Arm Resonance'},'FontSize', 16);
legend('PSL,VCO (in 26)', 'AOM,VCO (in 27)', 'EX,VCO (in 25)', 'COMM,VCO (in 28)', 'DIFF,VCO (in 29)', ...
        'Location', 'SouthEast')
ylabel('Mag [dB]', 'FontSize', 14);
axis([min(f) max(f) -300 20])
grid on
subplot(212);
semilogx(f, angle(a)*180/pi, ...
         f, angle(b)*180/pi, ...
         f, angle(c)*180/pi, ...
         f, angle(d)*180/pi, ...
         f, angle(e)*180/pi, 'LineWidth', 2)
ylabel('Phase [deg]', 'FontSize', 14);
xlabel('Freq [Hz]', 'FontSize', 14);
axis([min(f) max(f) 180*floor(min(angle(a)*180/pi)/180) 90*ceil(max(angle(a)*180/pi)/90)])
    set(gca,'YTick',[-1800:90:1800])
grid on
if save_figure == 1
    orient(hdl, 'landscape');
%    print('-dpdf', [save_figure_dir,num2str(hdl,'%.3d'), ' ',get(tt,'string'), '.pdf']);
    print('-dpdf', [save_figure_dir ,num2str(hdl,'%.3d'), '.pdf']);
end
