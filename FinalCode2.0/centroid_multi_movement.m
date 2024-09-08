close all; % Close all open figure windows
clear all; % Clear all variables from workspace
clc; % Clear command window

%% Initial setting
L = 1;         % Turbulence distance [m]
n = 1;         % Refractive index
dz = 0.2;      % Phase screen separation [m]
N = 2048;      % Sampling rate
D = 0.03;      % Phase screen size [m]
delta = D/N;   % Pixel spacing [m]
wavl = 1064e-9;  % Wavelength [m]
step_num = L/dz; % Number of phase screens
delta1 = delta;  % Set input pixel spacing
deltan = delta;  % Define output pixel spacing
k = n * 2 * pi / wavl;    % wavenumber
x = (-N/2:N/2-1) * delta; % spatial positions along x axis
y = x;
w0 = 3e-3; % Define radius parameter
[X, Y] = meshgrid(x, y); %form a 2D coordinate system over a plane
[theta, r] = cart2pol(X, Y); % Converts Cartesian coordinates into polar coordinates
CN = 4e-9;       % Set structure constant
R = 13;          % Define propagation distance

%% Simulation
% Generate flattened gaussian beam
M = 10;  % Beam order
E = zeros(N, N);
for m = 1:M
    E = E + (-1)^(m-1) * nchoosek(M, m) * exp(-m * (X.^2 + Y.^2) / w0^2);
end

% Apply lens
f1 = 1.5; % focus length
lens1_phase = lens_focus(D, N, wavl, f1);
EL = E .* lens1_phase; % After lens 1

% First 6m free-space propagation
zp1 = 6.1; % Propagation distance
[xn, yn, Ep1] = ang_spec_multi_prop_vac(EL, wavl, delta, delta, zp1);


% Number of simulations
num_simulations = 100;

% Preallocate arrays to store results
theta_x = zeros(step_num, num_simulations);
theta_y = zeros(step_num, num_simulations);
phz = zeros(N,N,step_num);


for i = 1:num_simulations
    % Generate phase screen
    E = Ep1;
    for m = 1:step_num
    phz(:,:,m) = vkolmg(D, dz, N, CN, wavl);
    E = E .* exp(1i * phz(:,:,m)); 
    I1 = abs(E).^2;
    I1 = I1 / max(I1(:));% Calculate intensity and normalize
    
    [~, ~, E] = ang_spec_multi_prop_vac(E, wavl, delta1, deltan, dz); 
    [theta_x(m,i), theta_y(m,i)] = cmove(I1,X,Y,dz); % calculate angular divergence
    end
end

%% Angular divergence of the beam centroid movement
figure(1);
hold on;  % Keep the plot

% Define a set of colors
colors = lines(step_num); 

% Initialize an empty array for legend names
legend_names = cell(1, step_num);

% Iterate through m in reverse, from step_num to 1
for m = step_num:-1:1
    plot(theta_x(m,:), theta_y(m,:), '-k', 'Color', colors(m, :)); 
    legend_names{step_num - m + 1} = sprintf('PS %d, z = %.1f', m, (6.1 + m*0.2));
    
    % Add labels and title
    xlabel('Horizontal deviation from mean [\murad]');
    ylabel('Vertical deviation from mean [\murad]');
    title('Beam Wander');
end

% After drawing all the lines, add the legend
legend(legend_names, 'Location', 'best');
hold off; 

%%  Temporal variation of the beam wander
figure(2);
hold on;  % Keep the plot

% Initialize another array for legend names
legend_names_2 = cell(1, step_num);

% Define the range of n (iterations)
n = 1:100;  

% Iterate through m in reverse, from step_num to 1
for m = step_num:-1:1
    % Plot the second graph, using different colors
    plot(n, theta_x(m,:), '-k', 'Color', colors(m, :)); 
    
    % Assign a name to each plot line
    legend_names_2{step_num - m + 1} = sprintf('PS %d, z = %.1f', m, (6.1 + m*0.2));
    
    % Add labels and title
    xlabel('n');
    ylabel('Horizontal deviation [\murad]');
    title('Horizontal Deviation Across n');
end

% After drawing all the lines, add the legend
legend(legend_names_2, 'Location', 'best');
hold off; 
