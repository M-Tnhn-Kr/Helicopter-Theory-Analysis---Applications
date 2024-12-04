%% Helicopter Theory HW1,Part2
clc; clear; close all;

%% Known Parameters (Hughes OH-6A "Cayuse")
rho = 1.225; % air density (kg/m^3)
omega = 485*(2*pi/60); % rotor axial speed (rpm converted to 1/s)  
R = 4.015; % rotor radius (m)
Vtip = omega * R; % tip speed (m/s)
A = pi*R^2; % disk area (m^2)
W_gr = 1089; % GTOW (kg)
g = 9.81; % gravitational acceleration (m/s^2)

%% Calculated Parameters
Ct_h = 1089*9.81/(rho*A*Vtip^2); % thrust coefficient for hover
lambda_h = sqrt(Ct_h/2); % induced speed at hover (will be used as initial value)

V_f = 0.01:1:100; % chosen forward flight speeds (m/s)
alpha_values = 0:2:8; % chosen disk tilt angles (degrees)
mu_values = V_f / Vtip; % advance ratios

%% Fixed Point Method Solution
figure; hold on;
for alpha = alpha_values
    % forward flight thrust coefficient
    Ct_f = W_gr * g * cosd(alpha) / (rho * A * Vtip^2);
    % inflow ratio (functions at last section)
    y = arrayfun(@(mu) fixed_point(Ct_f, mu, alpha, lambda_h), mu_values);
    % plotting
    plot(mu_values / lambda_h, y / lambda_h,"LineWidth",2.5);
end

% plot customization
legend("\alpha = 0","\alpha = 2","\alpha = 4","\alpha = 6","\alpha = 8","Location","northwest")
xlabel("Forward speed ratio (\mu / \lambda_h)");
ylabel("Inflow ratio (\lambda / \lambda_h)");
title("Variation of Inflow Ratio with Advance Ratio");

%% Newton Raphson Method
figure; hold on;
for alpha = alpha_values
    % Forward flight thrust coefficient
    Ct_f = W_gr * g * cosd(alpha) / (rho * A * Vtip^2);
    % Inflow ratio calculation using Newton-Raphson method
    y = arrayfun(@(mu) newton_raphson(Ct_f, mu, alpha, lambda_h), mu_values);
    % Plotting
    plot(mu_values / lambda_h, y / lambda_h, "LineWidth", 2.5);
end

% plot customization
legend("\alpha = 0","\alpha = 2","\alpha = 4","\alpha = 6","\alpha = 8","Location","northwest")
xlabel("Forward speed ratio (\mu / \lambda_h)");
ylabel("Inflow ratio (\lambda / \lambda_h)");
title("Variation of Inflow Ratio with Advance Ratio");

%% Used Functions

% 1. Fixed Point Iteration Function
function lambda = fixed_point(Ct_f, mu, alpha, lambda_h)
    lambda = lambda_h; % initial guess
    while true
        lambda_new = mu * tand(alpha) + Ct_f / (2 * sqrt(mu^2 + lambda^2));
        if abs(lambda_new - lambda) < 1e-5, break; end
        lambda = lambda_new;
    end
end

% 2. Newton Raphson Method
function lambda = newton_raphson(Ct_f, mu, alpha, lambda_h)
    lambda = lambda_h; % Initial guess
    tol = 1e-10; % Convergence tolerance
    max_iter = 100; % Maximum number of iterations
    
    for iter = 1:max_iter
        % Governing equation for inflow ratio
        f = lambda - (mu * tand(alpha) + Ct_f / (2 * sqrt(mu^2 + lambda^2)));
        % Derivative of the governing equation
        df = 1 + (Ct_f * lambda) / (2 * (mu^2 + lambda^2)^(3/2));
        
        % Update step using Newton-Raphson
        lambda_new = lambda - f / df;
        
        % Check for convergence
        if abs(lambda_new - lambda) < tol
            lambda = lambda_new;
            return;
        end
        
        lambda = lambda_new; % Update value for next iteration
    end
    
    error("Newton-Raphson did not converge within %d iterations.", max_iter);
end

%% Mehmet Tunahan Kara / 110200028
% 