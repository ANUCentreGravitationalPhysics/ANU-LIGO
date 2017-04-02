%flatTOP = [0.4 0.4 1 1 2 2 3.4 3.4];

global w

if resprun == 0
    PDHresp;
end

clear ServoTM ServoPM ServoUIM ServoTOP S_TM S_PM S_UIM S_TOP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Feedback Servo Response, Open-loop

% Set overall gain in the closed loop response
if gLP ~= 0
  gLP = 1;%Should be 1
end


gALL = 2.75e3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TM Servo
% MEVANS - made zeros complex, added resgain, and increased UGF
fc = 80;

zTM = [ 0.4 0.4];%1
pTM = [ 2 1e3];
kTM = 5e4 * prod(pTM) / prod(zTM(2:end));

% Set the PM Servo ZPK
zzz = [-2*pi.*zTM];
ppp = [-2*pi.*pTM];
kkk = kTM;

setGainAtF = 1e-2;  % frequency in Hz
gainAtF = 1;     % gain at that frequency

% Set the TM Servo ZPK
zzz = [-2*pi.*zTM];
ppp = [-2*pi.*pTM];

ServoTest = zpk(zzz, ppp, 1);
kkk = gainAtF / abs(evalfr(ServoTest, 2i * pi * setGainAtF));

% Gain stage into the TM actuators
kdB = 0;  %135
TM_gain = 10^(kdB/20);

%preTM = 10^( 0 /20);

% Setting Simulink Blocks
ServoTM = zpk(zzz, ppp, kkk);
gTM = 0;
coTM = cutoff(fc, 4);

% Setting up complex TFs
%S_TM = ss2complex(ServoTM, w); 
% cfilter = ss2complex(coTM, w);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PM Servo
% MEVANS - made zeros complex, added resgain, and increased UGF
fc = 40;

prePM = 7e4;

zPM = [filtRes(1,3)' 1];
pPM = [1e-2 1e-2 1e3];%[ fpPM 1000 1000];

setGainAtF = 1e-2;  % frequency in Hz
gainAtF = 1;     % gain at that frequency

% Set the PM Servo ZPK
zzz = [-2*pi.*zPM];
ppp = [-2*pi.*pPM];

ServoTest = zpk(zzz, ppp, 1);
kkk = gainAtF / abs(evalfr(ServoTest, 2i * pi * setGainAtF));

% Gain stage into the PM actuators
kdB = 0;
PM_gain = 10^(kdB/20);

% Setting Simulink Blocks
ServoPM = zpk(zzz, ppp, kkk);

gPM = 0;
coPM = cutoff(fc, 4);

% Getting up complex TFs
%S_PM = ss2complex(ServoPM, w); 
%cfilter = ss2complex(coPM, w);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
