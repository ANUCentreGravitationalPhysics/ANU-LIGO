function [damper] = local5ss(gain)

% Total gain of around 55 dB at 1 Hz is OK


shadow	= 2000;		%response of shadow sensor	(V/m)





current	= 50;		% volts to current conversion	(V/A) - ohms!
coilmag	= 1;		% coil magnet coupling		(N/A)


% Differentiator up to 10 Hz
fm1 = zpk(0,-2*pi*10,1);

% Elliptic low pass at 6 Hz
[z,p,k] = ellip(2,2,11,2*pi*6,'s');
fm2 = zpk(z,p,k);

% Resonant gain at the lowest longitudinal frequency
fm3 = resgain(0.44,2,3);
fm4 = resgain(4.42,2,33);

dampy = fm1 * fm2 * fm3 * fm4;

[magg,phph] = bode(dampy,2*pi*1);

dampy = dampy * (500/magg) * gain;



damper = ss(dampy);


