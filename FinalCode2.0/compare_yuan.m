% This code uses an ideal flat-top beam to compare with the results in the 
% literature, thereby validating the feasibility of the model.
close all;
clear all;
clc;

%% Initial setting
L = 1;         % Turbulence distance [m]
n = 1.4;         % Refractive index
dz = 0.2;      % Phase screen separation [m]
N = 2048;      % Sampling rate
D = 0.05;      % Phase screen size [m]
delta = D/N;   % Pixel spacing [m]
wavl = 1064e-9;% Wavelength [m]
step_num = L/dz;% Number of phase screens
k = n * 2 * pi / wavl;
x = (-N/2:N/2-1) * delta;
y = x;
R = 5e-3;  % Radius of the circular aperture
[X, Y] = meshgrid(x, y);
[theta, r] = cart2pol(X, Y);
delta1 = delta;
deltan = delta; % Assuming the same grid spacing at the output
k= 4.7e-10;
CN0 = ones(1,5)*k;% structure constant


% Define figure
FigFinal(:,:) = figure;

%% Generate ideal flat-top beam (through a circular aperture)
A = (X.^2 + Y.^2) <= (R^2);
M = zeros(N,N);
M(A) = 1;


%% Apply lens
f1 = 50;
lens1_phase = lens_focus(D, N, wavl, f1);
EL = M .* exp(1i * lens1_phase);


%% First 6m free-space propagation
zp1 = 3; % Propagation distance
[xn, yn, Ep1] = ang_spec_multi_prop_vac(EL, wavl, delta, delta, zp1);


%% Propagation with multiple phase screens
figure; % Create one figure for all subplots
[xn, yn, E] = ang_spec_multi_prop_vac(Ep1, wavl, delta1, deltan, 0.1);

for m = 1:step_num
    phz(:,:,m) = vkolmg(D, dz, N, CN(m), wavl);
    E = E .* exp(1i * phz(:,:,m)); % Add a layer of random phase screen
    % Calculate intensity and normalize
    I1 = abs(E).^2;
    I1 = I1 / max(I1(:));

    % Plotting each result in a subplot
    subplot(2, ceil(step_num/2), m); % Adjust the subplot grid as needed
    plot_beam_D_C(xn(1,:), yn(:,1), I1, ...
        sprintf('Exit surface of Screen %d (z = %.1f)', m, 6.1+(m-1)*0.2));

    if m == 5
        [xn, yn, E] = ang_spec_multi_prop_vac(E, wavl, delta1, deltan, 0.1);
    else
        [xn, yn, E] = ang_spec_multi_prop_vac(E, wavl, delta1, deltan, dz);
    end
end

%% First turbulence region free-space propagation (SPF)
Es1 = Ep1;
[~, ~, Es1] = ang_spec_multi_prop_vac(Ep1, wavl, delta1, deltan, 0.5);

% Single phase screen
combined_ph = sum(phz, 3);
Es = Es1 .* exp(1i * combined_ph);
[xn, yn, Es] = ang_spec_multi_prop_vac(Es, wavl, delta1, deltan, 0.4);

Is = abs(Es).^2;
Is = Is / max(Is(:));

subplot(2,3,6);
plot_beam_D_C(xn(1,:), yn(:,1), Is, 'Exit surface of one Cumulative PS (z = 6.9)');

%% Final free-space propagation (MPF)
zp2 = 3; % Final propagation distance
[xn, yn, Emp] = ang_spec_multi_prop_vac(E, wavl, delta1, deltan, zp2);
Imp = abs(Emp).^2;
Imp = Imp / max(Imp(:));

figure(FigFinal(:,:)); subplot(1, 2, 1);
plot_beam_D_C(xn(1,:), yn(:,1), Imp, 'Final Intensity with Multiple PS (z = 13)');

%% Final free-space propagation (SPS)
Es2 = Es;
[~, ~, Es2] = ang_spec_multi_prop_vac(Es2, wavl, delta1, deltan, 3);
% Calculate intensity and normalize
Is2 = abs(Es2).^2;
Is2 = Is2 / max(Is2(:));

figure(FigFinal(:,:)); subplot(1, 2, 2);
plot_beam_D_C(xn(1,:), yn(:,1), Is2, 'Final Intensity with Single PS (z = 6.9)');



