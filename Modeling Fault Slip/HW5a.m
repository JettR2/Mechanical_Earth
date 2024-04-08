% Load the GPS displacement data from the Izmit_data.mat file
load('Izmit_data.mat');
north_south_km = Izmit_data(:, 1); % North-south distances in kilometers
westward_displacement_m = Izmit_data(:, 2); % Westward displacement in meters

% Define the fault slip and x1 values
s = 1; % Magnitude of Slip in meters
x1 = (north_south_km) * 1000; % NS slip values in kilometers, * 1000 to be in meters
d2 = 0; % Depth of the fault at the surface


% Define different fault depths (d1 values)
depths = [1000, 10000, 20000]; % 1km, 10km and 20km in meters for depth

% Initialize a figure
figure;

% Plot the predicted displacements for different depths
hold on;

for i = 1:length(depths)
    d1 = depths(i);
    u3_model = (-s/pi) * (atan(x1/d1)-atan(x1/d2)); %East west displacement in x3 direction, Positive is west
    
    scatter(x1, u3_model, 'Linewidth', 2, 'DisplayName', sprintf('Depth = %d km', depths(i) / 1000));
end

plot(x1, westward_displacement_m, 'b', 'LineWidth', 2,'DisplayName', sprintf('GPS Data')); % GPS data in blue

hold off;


% Label the axes and provide a legend
axis([-8.5*10^4,8*10^4,-2,2])
xlabel('x1 (meters)');
ylabel('u3 (meters)');
legend('Location','north');

% Title and grid (optional)
title('Vert Displacement vs. NS Distance for Different Fault Depths');
grid on;


%%
% Answer to Question 1a
% As the Depth, d1, Increases, the displacement at a particular distance
% from the fault decreases. So at 1km, 100000m from the fault, displacement
% is -0.5meters, or 0.5 meters to the East, when the fault reaches 20km,
% that displacemnt is -0.44 meters or 0.44 meters to the East.