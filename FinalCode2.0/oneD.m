close all;
clear all;
clc;

%% Initial setting
L = 1;         % Turbulence distance [m]
n = 1;         % Refractive index
dz = 0.2;      % Phase screen separation [m]
N = 2048;      % Sampling rate
D = 0.03;       % Phase screen size [m]
delta = D/N;   % Pixel spacing [m]
wavl = 1064e-9; % Wavelength [m]
step_num = L/dz; % Number of phase screens

% Structure constants for different groups
k= 9e-9;
CN0 = ones(1,5)*k;
CN(1,:) = CN0;
CN(2,:)  = [9e-9 4e-9 4e-9 4e-9 9e-9];
CN(3,:)  = [4e-9 4e-9 9e-9 4e-9 4e-9];


for j =1:1
    k = n*2*pi/wavl; % wavenumber
    x = (-N/2:N/2-1)*delta;
    y = x;
    w0 = 2e-3; % waist radius
    [X, Y] = meshgrid(x, y);
    delta1 = delta;
    deltan = delta; % Assuming the same grid spacing at the output

    % Define figure
    FigFinal(:,:,j) = figure; %MPS,SPS and NO turbulence plots
    FigNo(:,:, j) = figure;

    %% Generate flattened gaussian beam
    M = 10;  % Beam order
    E = zeros(N, N);
    for m = 1:M
        E = E + (-1)^(m-1) * nchoosek(M, m) * exp(-m * (X.^2 + Y.^2) / w0^2);
    end
    I = E .* conj(E);
    I = I / max(max(I));

    % figure;
    % plot_1D(x, I);

    % % plot initial beam
    % figure;
    % plot_beam_D_C(x, y, abs(E).^2, 'Initial Flattened Gaussian Beam');

    %% Apply lens
    f1 = 50;
    lens1_phase = lens_focus(D, N, wavl, f1);
    EL = E .* lens1_phase; % After lens 1
    IL = EL .* conj(EL);
    IL = IL / max(max(IL));

    % figure;
    % plot_beam_D_C(x, y, abs(EL).^2, 'Exit surface of Lens (z = 0)');



    %% First 6m free-space propagation
    zp1 = 6; % Propagation distance
    [xn, yn, Ep1] = ang_spec_multi_prop_vac(EL, wavl, delta, delta, zp1);
    Ip1 = Ep1 .* conj(Ep1);
    Ip1 = Ip1 / max(max(Ip1));


    % figure(FigNo(:,:,j)); subplot(2, 3, 1);
    % plot_beam_D_C(xn(1,:), yn(:,1), Ip1, 'Free-space propagation at z=6');



    %% Propagation with multiple phase screens
    figure; % Create one figure for all subplots
    [xn, yn, E] = ang_spec_multi_prop_vac(Ep1, wavl, delta1, deltan, 0.1);

    for m = 1:step_num
        phz(:,:,m) = vkolmg(D, dz, N, CN(j,m), wavl);
        E = E .* exp(1i * phz(:,:,m)); % Add a layer of random phase screen

        if m == 5
            [xn, yn, E] = ang_spec_multi_prop_vac(E, wavl, delta1, deltan, 0.1);
        else
            [xn, yn, E] = ang_spec_multi_prop_vac(E, wavl, delta1, deltan, dz);
        end

        % Calculate intensity and normalize
        I1 = abs(E).^2;
        I1 = I1 / max(I1(:));

        % Plotting each result in a subplot
        subplot(2, step_num+1, m); % Adjust the subplot grid as needed
        plot_beam_D_C_big(xn(1,:), yn(:,1), I1, sprintf('Exit surface of S%d (z = %.1f)', m, 6.1+(m-1)*0.2));

        % Plot 1-D
        subplot(2, step_num+1, m+step_num+1);
        plot_1D(xn(1,:), I1);

    end

    %% First turbulence region free-space propagation (SPF)

    [~, ~, Es1] = ang_spec_multi_prop_vac(Ep1, wavl, delta1, deltan, 0.5);

    % single phase screen
    combined_ph = sum(phz, 3);
    Es = Es1 .* exp(1i * combined_ph);
    [xn, yn, Es] = ang_spec_multi_prop_vac(Es, wavl, delta1, deltan, 0.4);

    Is = abs(Es).^2;
    Is = Is / max(Is(:));

    subplot(2, step_num+1, 6);
    plot_beam_D_C_big(xn(1,:), yn(:,1), Is, '     Exit surface of Cumulative PS (z = 6.9)');
    % Plot 1-D
    subplot(2, step_num+1, 12);
    plot_1D(xn(1,:), Is);




    %% Final free-space propagation (MPF)
    zp2 = 6; % Final propagation distance
    [xn, yn, Emp] = ang_spec_multi_prop_vac(E, wavl, delta1, deltan, zp2);
    Imp = abs(Emp).^2;
    Imp = Imp / max(Imp(:));

    figure(FigFinal(:,:,j)); subplot(2, 3, 1);
    plot_beam_D_C(xn(1,:), yn(:,1), Imp, 'Final Intensity with Multiple PS (z = 13)');
    subplot(2, 3, 4);
    plot_1D(xn(1,:), Imp);

    %% Final free-space propagation (SPS)

    Es2 = Es;
    [~, ~, Es2] = ang_spec_multi_prop_vac(Es2, wavl, delta1, deltan, 6.1);
    % Calculate intensity and normalize
    Is2 = abs(Es2).^2;
    Is2 = Is2 / max(Is2(:));

    figure(FigFinal(:,:,j));
    subplot(2, 3, 2);
    plot_beam_D_C(xn(1,:), yn(:,1), Is2, 'Final Intensity with Single PS (z = 13)');
    subplot(2, 3, 5);
    plot_1D(xn(1,:), Is2);



    %% NO turbulence 13m propagation
    E13 = Ep1;
    for i = 1:5
        [x13, y13, E13] = ang_spec_multi_prop_vac(E13, wavl, delta1, deltan, dz);
        I13 = E13 .* conj(E13);
        I13 = I13 / max(max(I13));

        % figure(FigNo(:,:,j)); subplot(2, 5, i);
        % plot_beam_D_C(x13(1,:), y13(:,1), I13, sprintf('Intensity without Turbulence at z = %.1f',6+i*dz));
        % % Plot 1-D
        % subplot(2, step_num, i+step_num);
        % plot_1D(x13(1,:), I13);


    end
    [x13, y13, E13] = ang_spec_multi_prop_vac(E13, wavl, delta1, deltan, 6);
    %[x13, y13, E13] = ang_spec_multi_prop_vac(E13, wavl, delta1, deltan, 7);

    I13 = E13 .* conj(E13);
    I13 = I13 / max(max(I13));

    % Plotting the results
    figure(FigFinal(:,:,j)); subplot(2, 3, 3);
    plot_beam_D_C(x13(1,:), y13(:,1), I13, 'Final Intensity without Turbulence (z = 13)');
    subplot(2, 3, 6);
    plot_1D(x13(1,:), I13);

end

% %% Plot Phase Screen
% for i=1:5
%     figure;
%     sgtitle('Phase Screen generated using Von Karman PSD');
%     subplot(121);imagesc(x, y,phz(:,:,i)),colorbar;
%     xlabel('x (m)'); ylabel('y (m)');
%     subplot(122);mesh(x, y,phz(:,:,2)),colorbar;
%     xlabel('x (m)'); ylabel('y (m)');zlabel('Intensity');
% end
