% flag to save the figures
% 1 = yes, 0 = not
save_figure = 0;

%%%%%%%%%%
% Plotting
figure(100)
hndl = loglog(f, abs(H1), ...
       f, abs(H2), ...
       f, abs(H5), ...
       f, abs(H6));
set(hndl, 'LineWidth', 2);
grid on;
title('Transfer Function - H = x/x')
xlabel('Frequency [Hz]');
ylabel('Magnitude [m/m]');
legend('GND -> TM', ...
       'PRN -> TM', ...
       'PRN -> RO', ...
       'GND -> RO', ...
       'Location', 'NorthEast');
axis([0.01 100 1e-8 1e0])

if save_figure == 1
    print -dpng sim/f100_Transfer_Function.png
end

% figure(110)
% subplot(2,1,1)
% hndl = semilogx(f, 20*log10(abs(S_PM)));
% set(hndl, 'LineWidth', 2);
% title('Servo Response - S_{PM}');
% ylabel('Magnitude [dB]');
% axis([1e-2 100 60 180]);
% grid on;
% legend(['Poles@',num2str(pp),'Hz, Zeros@', num2str(zz),'Hz, 4^{th}-order Bessel at ',num2str(fc_bessel),'Hz'], ...
%        'Location', 'NorthWest');
% subplot(2,1,2)
% hndl = semilogx(f, angle(S_PM).*180/pi);
% set(hndl, 'LineWidth', 2);
% ylabel('Phase [deg]');
% xlabel('Frequency [Hz]');
% axis([1e-2 100 -180 180]);
% grid on;
% 
% if save_figure == 1
%     print -dpng sim/f110_Servo_Response.png
% end

figure(120)
hndl = loglog(f, abs(H3), ...
              f, abs(H4), ...
              f, abs(X_PRN_TM));
axis([1e-2 100 1e-18 1e-5]);
grid on;
set(hndl, 'LineWidth', 2);
title('Transfer Function - H = [m/N]');
ylabel('Amplitude [m/N]');
xlabel('Frequency [Hz]');
legend('F_{TM} -> X_{TM}', ...
       'F_{PM} -> X_{TM}', ...
       'Location', 'SouthWest');

if save_figure == 1
    print -dpng sim/f120_Transfer_Function_.png
end

figure(140)
hndl = loglog(f, I_OSEM_PM, ...
      f, I_OSEM_PM_rms, '--b');
set(hndl, 'LineWidth', 2);
grid on;
%axis([1e-2 1e2 1e-8 1e0]);
title('Single OSEM Feedback Current');
ylabel('Current [A/rtHz]');
xlabel('Frequency [Hz]');
legend(['PM OSEMs - ',num2str(I_OSEM_PM_rms(1)*1000,3),' mA_{rms}'], ...
       'Location', 'SouthWest');

if save_figure == 1
    print -dpng sim/f140_Single_OSEM_fb_Current.png
end

figure(150)
hndl = loglog(f, abs(H_TM_Fpm), '-b');
set(hndl, 'LineWidth', 2);
grid on;
axis([1e-2 1e2 1e-10 1e-2]);
title('Total Force To Stage');
ylabel('Force [N/rtHz]');
xlabel('Frequency [Hz]');
legend('F_{PM}', ...
       'Location', 'SouthWest');

if save_figure == 1
    print -dpng sim/f150_Total_Force2Stage.png
end



figure(126)
subplot(2,1,1)
hndl = semilogx(f, 20*log10(abs(H7)), ...
                f, 20*log10(abs(H8)));
set(hndl, 'LineWidth', 2);
title('Transfer Functions');
ylabel('Magnitude [dB]');
% axis([1e-2 100 -100 150]);
grid on;
legend('X_{GND} -> F_{PM}', ...
       'X_{PRN} -> F_{PM}', ...
       'Location', 'NorthEast');
subplot(2,1,2)
hndl = semilogx(f, angle(H7).*180/pi, ...
                f, angle(H8).*180/pi);
set(hndl, 'LineWidth', 2);
ylabel('Phase [deg]');
xlabel('Frequency [Hz]');
axis([1e-2 100 -180 180]);
grid on;

if save_figure == 1
    print -dpng sim/f126_Servo_Response.png
end
