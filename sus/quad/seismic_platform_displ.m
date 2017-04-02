% Calculates the Seismic Platform Displacement Requirement
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seismic Platform Displacement Noise
% This is a noisy estimate (lots of people, say in 10 years time)
N_F = [0.01 0.03 0.1 0.15 0.2 0.3 0.5 1 10 30 100];
N_X = [5e-6 1e-6 1e-6 2e-6 5e-7 2e-7 8e-10 3e-11 6e-13 6e-14 3e-15];
sei_noisy_0612 = 10.^(interp1(N_F,log10(N_X),f,'cubic',-14));
X_sei_noisy = sei_noisy_0612;

% Seismic Platform Displacement Noise
% The is the current requirement
SEI_F = [0.01 0.03 0.1 0.2 0.5 1 10 30 100];
SEI_X = [3e-6 1e-6 2e-7 2e-7 8e-10 1e-11 3e-13 3e-14 1e-15];
sei_model_0612 = 10.^(interp1(SEI_F,log10(SEI_X),f,'cubic',-14));
X_sei_std = sei_model_0612;


clear N_F N_X SEI_F SEI_X;
