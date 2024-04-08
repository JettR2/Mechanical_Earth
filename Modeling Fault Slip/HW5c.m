% Load the GPS displacement data from the Izmit_data.mat file
load('Izmit_data.mat');
north_south_km = Izmit_data(:, 1); % North-south distances in kilometers
westward_displacement_m = Izmit_data(:, 2); % Westward displacement in meters

% Define the initial parameter guesses for s and d1
s_guess = 4; % Initial guess for slip
d1_guess = 5 * 1000; % Initial guess for depth (5 km)
d2 = 0;

% Perform parameter optimization to minimize differences
options = optimset('Display', 'iter'); % Enable display of optimization iterations
parameters = fminsearch(@(params) optimize_parameters(params, north_south_km, westward_displacement_m), [s_guess, d1_guess], options);


% Extract the optimized parameters
s_optimized = parameters(1);
d1_optimized = parameters(2);

% Print the values of s and d1 of the best fit
fprintf('Best fit parameters:\n');
fprintf('s (slip) = %.2f meters\n', s_optimized);
fprintf('d1 (depth) = %.2f meters\n', d1_optimized);

% Calculate the model predictions for the optimized parameters
x1 = north_south_km * 1000; % Convert north-south distances from kilometers to meters
u3_model = (-s_optimized/pi) * (atan(x1/d1_optimized)-atan(x1/d2));

% Create a plot of the model predictions and overlay with the GPS data
figure;

% Plot model predictions and GPS data
hold on;
plot(x1, westward_displacement_m, 'b', 'LineWidth', 2,'DisplayName', sprintf('GPS Data')); % GPS data in blue
scatter(x1, u3_model, 'r', 'LineWidth', 2,'DisplayName',sprintf('Best Fit Model')); % Model predictions in red
xlabel('x1 (meters)');
ylabel('u3 (meters)');
legend('Location','northwest');
title('Vertical Displacement vs. North-South Distance');
axis([-8.5*10^4,8*10^4,min([westward_displacement_m; u3_model]),max([westward_displacement_m; u3_model])]); % Adjust axis limits
grid on;
hold off;

function misfit = optimize_parameters(params, north_south_km, westward_displacement_m)
    s = params(1);
    d1 = params(2);
    d2 = 0; % Fixed value for d2
    
    % Calculate model predictions for the given parameters
    x1 = north_south_km * 1000; % Convert north-south distances from kilometers to meters
    u3_model = (-s/pi) * (atan(x1/d1)-atan(x1/d2));
    
    % Calculate the misfit as the sum of squared differences
    misfit = sum((u3_model - westward_displacement_m).^2);
end



%%
% We often oversimplify the fault geometry. We assume the fault is
% infinetly long, vertical strike slip with uniform slip over a finite
% depth in an elastic half space. This is a close enough approximation for
% the most part but it doesn't account for true fault geometry, disparities
% , anisotropy,etc. There may also be gps errors. 
