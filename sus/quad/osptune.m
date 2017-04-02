% QUADRESP makes some plots of the quad pendulum response

% Runs the quad model 4 different times to simulate 4 power levels


f_low = 0.1;
f_high = 200;

f = logspace(log10(f_low),log10(f_high),1000);
w = 2*pi*f;

% Set radiation pressure torque to zero
k_major=0;
k_minor=0;

k_ospring = 0;
damper = 1;
linquad;
quad_0W = ss_sys;

k_ospring = -1e5 * zpk([],-2*pi*[400],2*pi*400);
damper = 1;
linquad;
quad_1W = ss_sys;

k_ospring = -1e6 * zpk([],-2*pi*[400],2*pi*400);
damper = 1;
linquad;
quad_10W = ss_sys;

k_ospring = -1e7 * zpk([],-2*pi*[400],2*pi*400);
damper = 1;
linquad;
quad_125W = ss_sys;

bode(quad_0W(1,8),quad_1W(1,8),quad_10W(1,8),quad_125W(1,8),w)
title('Radiation Pressure -> Optic')
%axis([f_low f_high -160 -20])
legend('0 W','1 W','10 W','125 W')

save quadsys_rad quad_0W quad_1W quad_10W quad_125W

bob = quad_0W;

n = min(find(w > 2*pi*10));

[m0,p0] = bode(bob(1,18),w(n)); % TOP
[m1,p0] = bode(bob(1,10),w(n)); % UIM
[m2,p0] = bode(bob(1,9),w(n));  % PM
[m3,p0] = bode(bob(1,8),w(n));  % TM

meters_per_newton_10Hz = [m3;m2;m1;m0]


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









