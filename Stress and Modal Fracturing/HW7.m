% Question 1 A


% Parameters
sigma1 = 0;
sigma3 = 0;
P = 1; % Internal pressure
a = 1; % Half-length
b = 0.25; % Half-width
beta_range = 0:15:90; % Inclination angle range
eta_range = -90:1:90; % Position around the hole range

% Initialize the figure
figure;

% Plot the stress distributions for different inclinations
hold on;

for beta = beta_range
    % Calculate tangential stress
    tangential_stress = -P + ((sigma1 + sigma3 + 2*P) * 2*a*b ./ ((a^2 + b^2) - (a^2 - b^2) * cosd(2*eta_range))) - ...
        (((sigma1 - sigma3) * ((a + b)^2 * cosd(2*(beta - eta_range)) - (a^2 - b^2) * cosd(2*beta))) ./ ...
        ((a^2 + b^2) - (a^2 - b^2) * cosd(2*eta_range)));
    
    % Plot the tangential stress distribution
    plot(eta_range, tangential_stress, 'DisplayName', sprintf('Beta = %d Degrees', beta));
end

hold off;

% Label the axes and provide a legend
xlabel('Position around Hole (\eta)');
ylabel('Tangential Stress');
legend('show');

% Title and grid (optional)
title('Q1A) Tangential Stress Distribution for Different Inclinations');
grid on;

% Within the first set of parameters, there no remote stresses and an
% internal pressure of unit magnitude 1. This simplifies the tangential
% stress equation, because both numerators are a product of both remote
% stresses, so as they are set to zero, both numerators are zero. Since
% Beta, the inclination angle, is in the numerator, tangential stress does
% not vary with the inclination angle. 

%%
% Question 1 B

% Parameters
sigma1 = 1;
sigma3 = 0;
P = 0; % Internal pressure
a = 1; % Half-length
b = 0.25; % Half-width
beta_range = 0:15:90; % Inclination angle range
eta_range = -90:1:90; % Position around the hole range

% Initialize the figure
figure;

% First subplot: Tangential stress distribution for different inclinations
subplot(2,1,1); % 2 rows, 1 column, first subplot
hold on;

for beta = beta_range
    tangential_stress = -P + ((sigma1 + sigma3 + 2*P) * 2*a*b ./ ((a^2 + b^2) - (a^2 - b^2) * cosd(2*eta_range))) - ...
        (((sigma1 - sigma3) * ((a + b)^2 * cosd(2*(beta - eta_range)) - (a^2 - b^2) * cosd(2*beta))) ./ ...
        ((a^2 + b^2) - (a^2 - b^2) * cosd(2*eta_range)));
    
    plot(eta_range, tangential_stress, 'DisplayName', sprintf('Beta = %d Degrees', beta));
end

hold off;

xlabel('Position around Hole (\eta)');
ylabel('Tangential Stress');
legend('show');
title('Q1B) Tangential Stress Distribution for Different Inclinations');
grid on;

% Second subplot: Maximum tangential stress trend
subplot(2,1,2); % 2 rows, 1 column, second subplot

max_stresses = zeros(size(beta_range));
max_beta_values = zeros(size(beta_range));

for i = 1:length(beta_range)
    beta = beta_range(i);
    
    tangential_stress = -P + ((sigma1 + sigma3 + 2*P) * 2*a*b ./ ((a^2 + b^2) - (a^2 - b^2) * cosd(2*eta_range))) - ...
        (((sigma1 - sigma3) * ((a + b)^2 * cosd(2*(beta - eta_range)) - (a^2 - b^2) * cosd(2*beta))) ./ ...
        ((a^2 + b^2) - (a^2 - b^2) * cosd(2*eta_range)));
    
    [max_stress, max_index] = max(tangential_stress);
    max_eta_value = eta_range(max_index);
    
    max_stresses(i) = max_stress;
    max_beta_values(i) = beta;
end

plot(max_beta_values, max_stresses, 'o-', 'LineWidth', 2);
xlabel('Beta (Degrees)');
ylabel('Maximum Tangential Stress');
title('Max stress Values');
grid on;


%%

% Question 1 C

% Parameters
sigma1_c = 0;
sigma3_c = -1;
P_c = 0; % Internal pressure
a_c = 1; % Half-length
b_c = 0.25; % Half-width
beta_range_c = 0:15:90; % Inclination angle range
eta_range_c = -90:1:90; % Position around the hole range

% Initialize the figure
figure;

% First subplot: Tangential stress distribution for different inclinations (Question 1 C)
subplot(2,1,1); % 2 rows, 1 column, first subplot
hold on;

for beta = beta_range_c
    tangential_stress_c = -P_c + ((sigma1_c + sigma3_c + 2*P_c) * 2*a_c*b_c ./ ((a_c^2 + b_c^2) - (a_c^2 - b_c^2) * cosd(2*eta_range_c))) - ...
        (((sigma1_c - sigma3_c) * ((a_c + b_c)^2 * cosd(2*(beta - eta_range_c)) - (a_c^2 - b_c^2) * cosd(2*beta))) ./ ...
        ((a_c^2 + b_c^2) - (a_c^2 - b_c^2) * cosd(2*eta_range_c)));
    
    plot(eta_range_c, tangential_stress_c, 'DisplayName', sprintf('Beta = %d Degrees', beta));
end

hold off;

xlabel('Position around Hole (\eta)');
ylabel('Tangential Stress');
legend('show');
title('Q1C)Tangential Stress Distribution for Different Inclinations');
grid on;

% Second subplot: Maximum tangential stress trend (Question 1 C)
subplot(2,1,2); % 2 rows, 1 column, second subplot

max_stresses_c = zeros(size(beta_range_c));
max_beta_values_c = zeros(size(beta_range_c));

for i = 1:length(beta_range_c)
    beta = beta_range_c(i);
    
    tangential_stress_c = -P_c + ((sigma1_c + sigma3_c + 2*P_c) * 2*a_c*b_c ./ ((a_c^2 + b_c^2) - (a_c^2 - b_c^2) * cosd(2*eta_range_c))) - ...
        (((sigma1_c - sigma3_c) * ((a_c + b_c)^2 * cosd(2*(beta - eta_range_c)) - (a_c^2 - b_c^2) * cosd(2*beta))) ./ ...
        ((a_c^2 + b_c^2) - (a_c^2 - b_c^2) * cosd(2*eta_range_c)));
    
    [max_stress_c, max_index_c] = max(tangential_stress_c);
    max_eta_value_c = eta_range_c(max_index_c);
    
    max_stresses_c(i) = max_stress_c;
    max_beta_values_c(i) = beta;
end

plot(max_beta_values_c, max_stresses_c, 'o-', 'LineWidth', 2);
xlabel('Beta (Degrees)');
ylabel('Maximum Tangential Stress');
title('Max Stress Values');
grid on;



%%

% Question 1 D


% Parameters
sigma1 = 0;
sigma3 = -1;
P = 0; % Internal pressure
a = 1; % Half-length
b_range = [0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10]; % Half-width
beta = 60; % Inclination angle range
eta_range = -90:1:90; % Position around the hole range

% Initialize the figure
figure;

% Plot the stress distributions for different inclinations
hold on;

for b = b_range
    % Calculate tangential stress
    tangential_stress = -P + ((sigma1 + sigma3 + 2*P) * 2*a*b ./ ((a^2 + b^2) - (a^2 - b^2) * cosd(2*eta_range))) - ...
        (((sigma1 - sigma3) * ((a + b)^2 * cosd(2*(beta - eta_range)) - (a^2 - b^2) * cosd(2*beta))) ./ ...
        ((a^2 + b^2) - (a^2 - b^2) * cosd(2*eta_range)));
    
    % Plot the tangential stress distribution
    plot(eta_range, tangential_stress, 'DisplayName', sprintf('Width = %0.2f', b));
end

hold off;

% Label the axes and provide a legend
xlabel('Position around Hole (\eta)');
ylabel('Tangential Stress');
legend('show');

% Title and grid (optional)
title('Q1D) Tangential Stress Distribution for Different Cavity Widths');
grid on;

%%

% Question 2 A

% Define the polar angle range
theta = linspace(-pi, pi, 1000);

% Define the stress intensity factor (replace with your specific value)
K_I = 1;

% Calculate the stress components
sigma_rr = (K_I/(sqrt(2*pi))).*cos(theta/2).*(1+sin(theta/2).^2);
sigma_thetatheta = (K_I/(sqrt(2*pi))).*cos(theta/2).^3;
sigma_rtheta = (K_I/(sqrt(2*pi))).*sin(theta/2).*cos(theta/2).^2;

% Normalize by the stress at the crack tip
sigma_0 = sigma_rr(1); % Stress at the crack tip
sigma_rr_normalized = sigma_rr/sigma_0;
sigma_thetatheta_normalized = sigma_thetatheta/sigma_0;
sigma_rtheta_normalized = sigma_rtheta/sigma_0;

% Plot the normalized stress components
figure;
plot(theta, sigma_rr_normalized, 'r', 'LineWidth', 2, 'DisplayName', '\sigma_{rr}');
hold on;
plot(theta, sigma_thetatheta_normalized, 'g', 'LineWidth', 2, 'DisplayName', '\sigma_{\theta\theta}');
plot(theta, sigma_rtheta_normalized, 'b', 'LineWidth', 2, 'DisplayName', '\sigma_{r\theta}');
hold off;

% Label the axes and provide a legend
xlabel('\theta');
ylabel('Normalized Stress Components');
legend('Location', 'best');
title('Q2A) Mode I Stresses vs. \theta');
grid on;


%%

% Question 2 B

% Define the polar angle range
theta = linspace(-pi, pi, 1000);

% Define the stress intensity factor (replace with your specific value)
K_II = 1;

% Calculate the stress components and make Left lateral motion as
% positive values
sigma_rr = (-K_II/(sqrt(2*pi))).*sin(theta/2).*(1-3*sin(theta/2).^2);
sigma_thetatheta = (-K_II/(sqrt(2*pi))).*(-3*sin(theta/2)).*cos(theta/2).^2;
sigma_rtheta = (-K_II/(sqrt(2*pi))).*(1-3*sin(theta/2).^2).*cos(theta/2);

% Normalize by the stress at the crack tip
sigma_0 = sigma_rr(1); % Stress at the crack tip
sigma_rr_normalized = sigma_rr/sigma_0;
sigma_thetatheta_normalized = sigma_thetatheta/sigma_0;
sigma_rtheta_normalized = sigma_rtheta/sigma_0;

% Plot the normalized stress components
figure;
plot(theta, sigma_rr_normalized, 'r', 'LineWidth', 2, 'DisplayName', '\sigma_{rr}');
hold on;
plot(theta, sigma_thetatheta_normalized, 'g', 'LineWidth', 2, 'DisplayName', '\sigma_{\theta\theta}');
plot(theta, sigma_rtheta_normalized, 'b', 'LineWidth', 2, 'DisplayName', '\sigma_{r\theta}');
hold off;

% Label the axes and provide a legend
xlabel('\theta');
ylabel('Normalized Stress Components');
legend('Location', 'best');
title('Q2B) Mode II Stresses (Left Lateral) vs. \theta');
grid on;
