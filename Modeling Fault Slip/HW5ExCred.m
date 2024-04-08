% Parameters
mu = 1; % Shear modulus (arbitrary units)
s = 1; % Slip (arbitrary units)
d = 1; % Distance to the dislocations (arbitrary units)

% Nondimensional distance range
x2_range = linspace(-3*d, 3*d, 1000);

% Calculate sigma_13
sigma13 = (s * mu * d) ./ (pi * (x2_range - d) .* (x2_range + d));

% Calculate nondimensional stress
nondimensional_stress = (2 * pi * d * sigma13) / (mu * s);

% Initialize the figure
figure;

% Plot the nondimensional stress
plot(x2_range/d, nondimensional_stress, 'LineWidth', 2);

% Label the axes
xlabel('Nondimensional Distance (x_2 / d)');
ylabel('Nondimensional Stress (stress = 2*\pi*d*\sigma_{13} / \mu * s)');

% Title and grid (optional)
title('Nondimensional Stress vs. Nondimensional Distance');

% Display the plot
grid on;
