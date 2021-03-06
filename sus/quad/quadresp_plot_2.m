% flag to save the figures
% 1 = yes, 0 = not
save_figure = 0;

%%%%%%%%%%
% Plotting
figure(100)
subplot(2,1,1)
hndl = semilogx(f, 20*log10(abs(H1)), ...
                f, 20*log10(abs(H2)), ...
                f, 20*log10(abs(H5)), ...
                f, 20*log10(abs(H6)), ...
                f, 20*log10(abs(H9)));
set(hndl, 'LineWidth', 2);
grid on;
title('PM Transfer Function - H = x/x')
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
legend('GND -> TM', ...
       'PRN -> TM', ...
       'PRN -> RO,pm', ...
       'GND -> RO,pm', ...
       'H_{OL}= RO,pm -> TM', ...
       'Location', 'SouthWest');
axis([0.01 100 -200 1e2])
subplot(2,1,2)
hndl = semilogx(f, angle(H1).*180/pi, ...
                f, angle(H2).*180/pi, ...
                f, angle(H5).*180/pi, ...
                f, angle(H6).*180/pi, ...
                f, angle(H9).*180/pi);
set(hndl, 'LineWidth', 2);
ylabel('Phase [deg]');
xlabel('Frequency [Hz]');
axis([1e-2 100 -180 180]);
grid on;

if save_figure == 1
    print -dpng sim/f100_Transfer_Function.png
end

figure(101)
subplot(2,1,1)
hndl = semilogx(f, 20*log10(abs(H11)), ...
                f, 20*log10(abs(H12)), ...
                f, 20*log10(abs(H15)), ...
                f, 20*log10(abs(H16)), ...
                f, 20*log10(abs(H19)));
set(hndl, 'LineWidth', 2);
grid on;
title('UIM Transfer Function - H = x/x')
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
legend('GND -> PM', ...
       'PRN -> PM', ...
       'PRN -> RO,uim', ...
       'GND -> RO,uim', ...
       'H_{OL}= RO,uim -> PM', ...
       'Location', 'SouthWest');
axis([0.01 100 -200 1e2])


figure(110)
hndl = loglog(f, abs(H3), ...
              f, abs(H4));
%axis([1e-2 100 1e-18 1e-5]);
grid on;
set(hndl, 'LineWidth', 2);
title('Transfer Function - H = [m/N]');
ylabel('Amplitude [m/N]');
xlabel('Frequency [Hz]');
legend('F_{TM} -> X_{TM}', ...
       'F_{PM} -> X_{TM}', ...
       'Location', 'SouthWest');

if save_figure == 1
    print -dpng sim/f110_Transfer_Function_.png
end

figure(120)
hndl = loglog(f, (abs(H8)), ...
                f, (abs(H10)), ...
                f, abs(H18), ...
                f, abs(H20));
set(hndl, 'LineWidth', 2);
title('Displacement To Force Response');
ylabel('Amplitude [N/m]');
% axis([1e-2 100 -100 150]);
grid on;
legend('X_{PRN} -> F_{PM}', ...
       'X_{TM} -> F_{PM}', ...
       'X_{PRN} -> F_{UIM}', ...
       'X_{PM} -> F_{UIM}', ...
       'Location', 'SouthWest');

if save_figure == 1
    print -dpng sim/f120_Displ_To_Force.png
end

figure(130)
hndl = loglog(f, FF_PM, '-b', ...
              f, FF_UIM, '-g', ...
              f, FF_TM, '-m', ...
              f, FF_PM_rms, '--b', ...
              f, FF_UIM_rms, '--g', ...
              f, FF_TM_rms, '--m');
set(hndl, 'LineWidth', 2);
grid on;
axis([1e-2 1e2 1e-10 1e-2]);
title('Total Force To Stage');
ylabel('Force [N/rtHz]');
xlabel('Frequency [Hz]');
legend(['F_{PM} - ',num2str(FF_PM_rms(1)*1000,3),' mN_{rms}'], ...
       ['F_{UIM} - ',num2str(FF_UIM_rms(1)*1000,3),' mN_{rms}'], ...
       ['F_{TM} - ',num2str(FF_TM_rms(1)*1e6,3),' \muN_{rms}'], ...
       'Location', 'SouthWest');

if save_figure == 1
    print -dpng sim/f130_Total_Force2Stage.png
end

figure(140)
hndl = loglog(f, I_OSEM_PM, '-b', ...
              f, I_OSEM_UIM, '--g', ...
              f, I_OSEM_PM_rms, '--b', ...
              f, I_OSEM_UIM_rms, '--g');
set(hndl, 'LineWidth', 2);
grid on;
%axis([1e-2 1e2 1e-8 1e0]);
title('Single OSEM Feedback Current');
ylabel('Current [A/rtHz]');
xlabel('Frequency [Hz]');
legend(['PM OSEMs - ',num2str(I_OSEM_PM_rms(1)*1000,3),' mA_{rms}'], ...
       ['UIM OSEMs - ',num2str(I_OSEM_UIM_rms(1)*1000,3),' mA_{rms}'], ...
       'Location', 'SouthWest');

if save_figure == 1
    print -dpng sim/f140_Single_OSEM_fb_Current.png
end

figure(150)
hndl = loglog(f, X_TM_CL, '-r', ...
              f, X_GND_TM, '-m', ...
              f, X_PRN_TM, '-g', ...
              f, X_TM_CL_rms, '--r', ...
              f, X_GND_TM_rms, '--m', ...
              f, X_PRN_TM_rms, '--g');
set(hndl, 'LineWidth', 2);
grid on;
axis([1e-2 1e2 1e-15 1e-6]);
title('Diplacement Noise');
ylabel('Displacement [m/rtHz]');
xlabel('Frequency [Hz]');
legend('X_{TM} - Closeloop', ...
       'X_{TM} - Openloop', ...
       'X_{PRN-TM} - Closeloop', ...
       'Location', 'SouthWest');

if save_figure == 1
    print -dpng sim/f150_Displacement_Noise.png
end


figure(200)
subplot(2,1,1)
hndl = semilogx(f, 20*log10(abs(S_PM)), ...
                f, 20*log10(abs(S_UIM)));
set(hndl, 'LineWidth', 2);
title('Servo Response ');
ylabel('Magnitude [dB]');
axis([1e-2 100 60 180]);
grid on;
legend(['S_{PM}-P@',num2str(pPM),'Hz, Z@', num2str(zPM),'Hz, 4^{th}-order Bessel at ',num2str(fc_bessel),'Hz'], ...
       ['S_{UIM}-P@',num2str(pUIM),'Hz, Z@', num2str(zUIM),'Hz, 2^{th}-order Bessel at ',num2str(fc_bessel),'Hz'], ...
       'Location', 'NorthWest');
subplot(2,1,2)
hndl = semilogx(f, angle(S_PM).*180/pi, ...
                f, angle(S_UIM).*180/pi);
set(hndl, 'LineWidth', 2);
ylabel('Phase [deg]');
xlabel('Frequency [Hz]');
axis([1e-2 100 -180 180]);
grid on;

% if save_figure == 1
%     print -dpng sim/f110_Servo_Response.png
% end






