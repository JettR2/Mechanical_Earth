% Load the GPS displacement data from the Izmit_data.mat file
load('Izmit_data.mat');
north_south_km = Izmit_data(:, 1); % North-south distances in kilometers
westward_displacement_m = Izmit_data(:, 2); % Westward displacement in meters

% Define the fault slip, Fault Depth and x1 values
s = 1; % Slip in meters
d1 = 10000; % Fault depth in meters (10 km)
x1 = (north_south_km) * 1000; % NS displacemnt in meters
d2 = 0;

% Define different fault slips (s values)
slips = [1, 4, 7]; % Slip values in meters

% Initialize a figure
figure;

% Plot the predicted displacements for different slips
hold on;

for i = 1:length(slips)
    s = slips(i);
    u3_model = (-s/pi) * (atan(x1/d1)-atan(x1/d2));
    
    scatter(x1, u3_model, 'LineWidth', 2, 'DisplayName', sprintf('Slip = %d m', slips(i)));
end

plot(x1, westward_displacement_m, 'b', 'LineWidth', 2,'DisplayName', sprintf('GPS Data')); % GPS data in blue
hold off;

% Label the axes and provide a legend
axis([-8.5*10^4,8*10^4,-3.5,3.5])
xlabel('x1 (meters)');
ylabel('u3 (meters)');
legend('Location', 'NorthWest');

% Title and grid (optional)
title('Vert Displacement vs. NS Distance for Different Slips');
grid on;


%%
% Increasing the magnitude of the slip also increases the rate of 
% displacemnt until a distance of 100000m. Then all of the different scenarios
% reach a near 0 slope, but at different magnitudes of u3.