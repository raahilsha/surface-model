clear all;
close all; 

%% Define surface
N1 = 150;
N2 = 150;

z = zeros(N1, N2);
x = linspace(-5,5, N1);
y = linspace(-5,5, N2);

for i = 1:N1
    for j = 1:N2
        z(i,j) = -1 * sqrt(64 - x(i) ^ 2 - y(j) ^ 2);
    end
end

figure(1);
plots = surf(x, y, z');
set(plots, 'edgecolor', 'none');

% input ball's position
gx = 2.5;
gy = 2.5;
gz = -1 * sqrt(64 - gx ^ 2 - gy ^ 2);

m = 0.01; % kg
r = 0.01; % m
g = 9.8; % m/s^2

dt = 0.005; % s

wx = 70;
wy = 0;
time = 0;

maxN = 5000;
xtravel = linspace(-4,4, maxN);
ytravel = linspace(-4,4, maxN);
ztravel = linspace(-4,4, maxN);
timetravel = linspace(-4,4, maxN);

for N = 1:maxN
    time = time + dt;
    
    theta_x = atan(-1 * sqrt(64 - gx ^ 2 - gy ^ 2) / gx);
    theta_y = atan(-1 * sqrt(64 - gx ^ 2 - gy ^ 2) / gy);
    
    ax = 5 * g * sin(theta_x) / (2 * r);
    ay = 5 * g * sin(theta_y) / (2 * r);
    
    wx = wx + ax * dt;
    wy = wy + ay * dt;
    vx = r * wx;
    vy = r * wy;
    
    delta_x = vx * dt;
    delta_y = vy * dt;
    
    gx = gx + delta_x;
    gy = gy + delta_y;
    gz = -1 * sqrt(64 - gx ^ 2 - gy ^ 2);
    
    xtravel(N) = gx;
    ytravel(N) = gy;
    ztravel(N) = gz;
    timetravel(N) = time;
end

hold on
plot3(xtravel, ytravel, ztravel, 'Color', 'magenta', 'LineWidth', 3);