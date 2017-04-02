function [damper] = local4ss(gain)
% local dev from Stuart Killbourn by KS 6/98 for triplep LIGO2
% modified again for 30 kg triple pendulum 6/99 KAS.
%
% scale factor
% 5 mm range  = <10V so >2000 V/m

shadow	= 2000;		%response of shadow sensor	(V/m)
[ag,bg,cg,dg] = zp2ss([],[],shadow*gain);	%(V/m)

fcut	= 9;
dcGain	= 1;

[alp,blp,clp,dlp] = lowpass(fcut,dcGain);


fpeak	= 15;
qp	= 5;
fnotch	= 23;
qn	= 20;
dcGain	= 1.5;

% [an3,bn3,cn3,dn3] = sculte(fpeak,qp,fnotch,qn,dcGain);

% [an,bn,cn,dn] = series(alp,blp,clp,dlp,an3,bn3,cn3,dn3);

current	= 50;		% volts to current conversion	(V/A) - ohms!
coilmag	= 5;		% coil magnet coupling		(N/A)

[ad,bd,cd,dd] = zp2ss([],[],(coilmag/current));	%(N/V)

[ahp,bhp,chp,dhp] = highpass(0.7,1);

% stuarts with modification for more hf phase margin
% transdif is just a little lead filter
[ad1,bd1,cd1,dd1] = transdif(0.35,0.7,1);
[ad2,bd2,cd2,dd2] = transdif(1,9,1);

[al,bl,cl,dl] = series(ahp,bhp,chp,dhp,ad1,bd1,cd1,dd1);
[al,bl,cl,dl] = series(al,bl,cl,dl,ad2,bd2,cd2,dd2);
[al,bl,cl,dl] = series(al,bl,cl,dl,alp,blp,clp,dlp);
[al,bl,cl,dl] = series(al,bl,cl,dl,ad,bd,cd,dd);
[al,bl,cl,dl] = series(al,bl,cl,dl,ag,bg,cg,dg);

damper = ss(al,bl,cl,dl);


