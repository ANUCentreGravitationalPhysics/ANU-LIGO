% ANGRESP makes some plots of the quad pendulum response

% Assumes 'generate_simulink.m' has been run to generate the state-space
% model, 'quad_sys', of the quad pendulum

k_ospring = 0;

k_major = 0;
k_minor = 0;


f_low = 0.01;
f_high = 100;

f = logspace(log10(f_low),log10(f_high),1000);
w = 2*pi*f;


damper = 1;
linquad;
eddie = ss_sys;
damper = 4;
linquad;
free = ss_sys;


figure(6)
%set(hndl, 'LineWidth',2);
bodemag(eddie(5,5),eddie(1,5), free(1,5),w)
title('Quad Torque -> Optic')
legend('Ground Pitch -> Optic Pitch (ECD)', ...
        'Ground Pitch -> Optic Long (ECD)', ...
        'Ground Pitch -> Optic Long (free)', ...
    'Location','SouthWest')
axis([f_low f_high -140 20])
grid on
[mag_pit2optLong_free(1,:), phs_pit2optLong_free(1,:)] = bode(free(1,5),w);
[mag_pit2optLong_eddie(1,:), phs_pit2optLong_eddie(1,:)] = bode(eddie(1,5),w);
save quad_pit2long_free f mag_pit2optLong_free phs_pit2optLong_free;
save quad_pit2long_eddie f mag_pit2optLong_eddie phs_pit2optLong_eddie;
print -dpng torque2opt.png

figure(4)
bodemag(eddie(5,11),eddie(5,12),eddie(5,13),eddie(5,20),eddie(5,5),w)
title('Torque -> Optic Pitch')
legend('TM','PM','UIM','Top','Ground','Location','SouthWest')
axis([f_low f_high -100 20])

[m0,p0] = bode(eddie(5,11),w(1));
[m1,p0] = bode(eddie(5,12),w(1));
[m2,p0] = bode(eddie(5,13),w(1));
[m3,p0] = bode(eddie(5,20),w(1));

rad_per_torque_dc_pit = [m3;m2;m1;m0]

figure(5)
bodemag(eddie(6,15),eddie(6,16),eddie(6,17),eddie(6,19),w)
title('Torque -> Optic Yaw')
legend('TM','PM','UIM','Top','Location','SouthWest')
axis([f_low f_high -100 20])

[m0,p0] = bode(eddie(6,15),w(1));
[m1,p0] = bode(eddie(6,16),w(1));
[m2,p0] = bode(eddie(6,17),w(1));
[m3,p0] = bode(eddie(6,19),w(1));

rad_per_torque_dc_yaw = [m3;m2;m1;m0]

[mag(1,:), phs(1,:)] = bode(eddie(5,1),w);
H_long2pit = mag .* exp(i.*phs.*pi/180);

[mag(1,:), phs(1,:)] = bode(eddie(5,5),w);
H_pit2pit = mag .* exp(i.*phs.*pi/180);

[mag(1,:), phs(1,:)] = bode(eddie(6,6),w);
H_yaw2yaw = mag .* exp(i.*phs.*pi/180);

seismic_platform_displ;
TMpit1 = X_sei_std .* H_long2pit;
TMpit2 = X_sei_std .* H_pit2pit ./ 10;
TMyaw = X_sei_std .* H_yaw2yaw ./ 10;

figure(1)
loglog(f, abs(TMpit1), f, abs(TMpit2), f, abs(TMyaw));
ylabel('rad/rtHz');
xlabel('Frequency [Hz]');

figure(7)
bodemag(eddie(5,1), eddie(5,5), eddie(6,6), w);
title('Transfer Functions - GND->TM');
legend('long -> pitch', 'pitch -> pitch', 'yaw -> yaw');
%bodemag(geod(1,1),eddie(1,1),damp3(1,1),free_quad(1,1),w)
%save eddie.mat eddie;
