% SIDLES makes some plots of the quad pendulum response
%
% This includes the optical torques from radiation pressure
% Runs the quad model 4 different times to simulate 4 power levels


f_low = 0.01;
f_high = 40;

f = logspace(log10(f_low),log10(f_high),2000);
w = 2*pi*f;

k_ospring = 0;

k_major = 0;
k_minor = 0;
damper = 3;
linquad;
quad_0W = ss_sys;

k_major = -300 * zpk([],-2*pi*300,2*pi*300);
k_minor = 1;
damper = 3;
linquad;
quad_major = ss_sys;

k_major = 11;
k_minor = 1 * zpk([],-2*pi*300,2*pi*300);
damper = 3;
linquad;
quad_minor = ss_sys;

g = 1000;

mybode(quad_0W(6,15)*g,quad_major(6,15)*g,quad_minor(6,15)*g,w)
title('Torque -> Yaw Response')
%axis([f_low f_high -160 -20])
legend('0 W','major','minor')
orient tall
print -dpdf yawspring_minus.pdf
%save quadsys_rad quad_0W quad_1W quad_10W quad_125W

return

subplot(221)
bodemag(eddie(1,1),odor(1,1),w)
title('Ground -> Optic')
%legend('Eddy','GEO','Bobo','Free')
axis([f_low f_high -260 50])


subplot(223)
bodemag(eddie(1,7),odor(1,7),w)
title('Sensor Noise -> Optic')
axis([f_low f_high -260 20])

subplot(222)
bodemag(eddie(1,8),odor(1,8),w)
title('Radiation Pressure -> Optic')
axis([f_low f_high -160 -20])









