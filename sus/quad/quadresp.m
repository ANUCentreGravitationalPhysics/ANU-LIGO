%% QUADRESP makes some plots of the quad pendulum response

% Assumes 'generate_simulink.m' has been run to generate the state-space
% model, 'quad_sys', of the quad pendulum



f_low = 0.01;
f_high = 30;

f = logspace(log10(f_low),log10(f_high),250);
w = 2*pi*f;

% Angular radiation pressure torque coefficients
k_major = 0;
k_minor = 0;

%%
damper = 1;
linquad;
quadEddie = ss_sys;

%%
damper = 2;
linquad;
quadGEO = ss_sys;

%%
damper = 3;
linquad;
bobo = ss_sys;
quadFancy = ss_sys;
%save quadsys quad_sys

%%
damper = 4;
linquad;
quadNoDamping = ss_sys;
[mag_qd_free(1,:), phs_qd_free(1,:)] = bode(quadNoDamping(1,1),w);
%save quad_free f mag_qd_free phs_qd_free;

%save quadDamping quadEddie quadGEO quadFancy quadNoDamping;

%%
figure(1)
%subplot(221)
bodemag(quadEddie(1,1),quadGEO(1,1),quadFancy(1,1),quadNoDamping(1,1),w)
%bodemag(eddie(1,7),geod(1,7),bobo(1,7),free_quad(1,7),w)
grid on;
title('Ground -> Optic')
legend('Eddy','GEO','RG','Free')
axis([0.1 60 -160 50])
print -dpng g2opt.png
print -depsc g2opt.eps

%subplot(222)
figure(2)
bodemag(quadEddie(1,7),quadGEO(1,7),quadFancy(1,7),quadNoDamping(1,7),w)
title('Sensor Noise -> Optic')
legend('Eddy','GEO','RG','Free')
axis([0.1 20 -170 20])
print -dpng s2opt.png

%subplot(223)
figure(3)
bodemag(quadEddie(1,8),quadGEO(1,8),quadFancy(1,8),quadNoDamping(1,8),w)
title('Radiation Pressure -> Optic')
legend('Eddy','GEO','RG','Free')
axis([0.1 20 -120 0])
print -dpng r2opt.png

%%
figure(4)
bodemag(quadEddie(5,1),quadGEO(5,1),quadFancy(5,1),quadNoDamping(5,1),w)
title('GND -> Pit TM')
legend('Eddy','GEO','RG','Free')
axis([0.1 20 -60 40])

%bodemag(geod(1,1),eddie(1,1),damp3(1,1),free_quad(1,1),w)
