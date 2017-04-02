% PDHservo_All_2010_07_19_15_51_34
%
% TM and PM feedback servo settings
%
% J. Miller
%
%flatTOP = [0.4 0.4 1 1 2 2 3.4 3.4];

global w

% if resprun == 0
%     PDHresp;
% end
% 
% clear ServoTM ServoPM ServoUIM ServoTOP S_TM S_PM S_UIM S_TOP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Feedback Servo Response, Open-loop

% Set overall gain in the closed loop response
if gLP ~= 0
  gLP = 1;%Should be 1
end


gALL = 1.75e4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TM Servo
% MEVANS - made zeros complex, added resgain, and increased UGF
fc = 400;

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
fc = 75;

prePM = 1e5;

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
% UIM Servo
fc = 3;

%gALL = 1000;
preUIM = 3e2; % gain prior the UIM servo

fpUIM = 0.01;

zUIM = [2  2 2];%[fpPM 2 2];
pUIM = [1e-2 1 1e3];%[fpUIM 1000 1000];

setGainAtF = 1e-2;  % frequency in Hz
gainAtF = 1;     % gain at that frequency

% Set the Simulink Blocks
zzz = [-2*pi.*zUIM];
ppp = [-2*pi.*pUIM];

ServoTest = zpk(zzz, ppp, 1);
kkk = gainAtF / abs(evalfr(ServoTest, 2i * pi * setGainAtF));

% Gain stage into the UIM actuators
kdB = 0;
UIM_gain = 10^(kdB/20);

% Setting Simulink Blocks
ServoUIM = zpk(zzz, ppp, kkk);

gUIM = 0;
coUIM = cutoff(fc, 4);

% Getting up complex TFs
%S_UIM = ss2complex(ServoUIM, w);
%cfilter = ss2complex(coUIM, w);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOP Servo
fc = 1;
fpTOP = 0.001;

%zTOP = [0.45*r1];
zTOP = [ 3.4 3.4];%2 2 [fpUIM 3.4 3.4];    % [Hz]
pTOP = [1e-2 1e3];%[fpTOP 2000 2000]; % [Hz]

setGainAtF = 1e-2;  % frequency in Hz
gainAtF = 1;     % gain at that frequency

% Set the Simulink Blocks
zzz = [-2*pi.*zTOP];
ppp = [-2*pi.*pTOP];

ServoTest = zpk(zzz, ppp, 1);
kkk = gainAtF / abs(evalfr(ServoTest, 2i * pi * setGainAtF));

% Gain stage into the TOP actuators
kdB = 0;
TOP_gain = 10^(kdB/20); % becomes gTOP

preTOP = 1e2;

% Setting Simulink Blocks
ServoTOP = zpk(zzz, ppp, kkk);
gTOP = 0;
coTOP = cutoff(fc, 4);

% Getting up complex TFs
%S_TOP = ss2complex(ServoTOP, w);
%cfilter = ss2complex(coTOP, w);

