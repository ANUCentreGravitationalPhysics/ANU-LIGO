% 
% r = filtRes(f0, Q)
%
% f0 and Q are the resonance frequency and quaity factor
% r is the corresponding pair of roots (i.e., zeros or poles)

function r = filtRes(f0, Q)

  r0 = (f0 ./ (2 * Q)) .* (1 + sqrt(1 - 4 * Q.^2));
  r = [r0; conj(r0)];

% phi = atan(sqrt(4 * Q.^2 - 1));
% r = f0 * [exp(i * phi); exp(-i * phi)];
