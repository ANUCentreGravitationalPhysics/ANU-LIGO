% local damping functions for the quad pendulum model
% P Fritschel, 6 aug 2001
% standard local control added as option K Strain 28th Aug 2001
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                   %
%    eddy current damping model, proportional to    %
%    velocity up to ~150 Hz                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% gain included to represent appropriate eddy current coupling
% but why does this cutoff at ~150 Hz ?
veldamp = zpk([0],-1000,30000);  

if damper == 1
  damping  = veldamp;   % velocity damping (eddy current)
  
  ldamp = damping;
  pdamp = damping;
  vdamp = damping;
  ydamp = damping;
  tdamp = damping;
  rdamp = damping;


elseif damper == 2
  damping  = local4ss(1);   % standard adapted GEO active damping, normalized to 'unit' gain

  ldamp = damping;
  pdamp = damping;
  vdamp = damping;
  ydamp = damping;
  tdamp = damping;
  rdamp = damping;

elseif damper == 3
  %damping = local5ss(1);

  ldamp = local5ss(1);
  pdamp = local4ss(1);
  vdamp = local4ss(1);
  ydamp = local5ss(1);
  tdamp = local4ss(1);
  rdamp = local4ss(1);
  
elseif damper == 4

  damping = 1;   % no damping
  
  ldamp = damping;
  pdamp = damping;
  vdamp = damping;
  ydamp = damping;
  tdamp = damping;
  rdamp = damping;
  
end
% start simple: all damping paths have same damping function


