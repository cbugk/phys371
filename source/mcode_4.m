%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Celil Bugra Karacan - 21401700 - PHYS371 Bilkent
%Homework-IV Question-I Due: 5th Nov 2018
clear; clc;

syms z; % for analytical solutions

k = 1; % N/m
m = 1.2; % kg
g = 9.81; %m/s^2

N = 1000; % Number of steps
T = 50; % Total time enlapsed (s)

y0 = -10; % initial conditions
v0 = 0;

yAnalytic = y0*cos(sqrt(k/m)*z); % analytical solutions defined symbolically
vAnalytic = - sqrt(k/m)*y0*sin(sqrt(k/m)*z);

[y, v] = NumericalFun(y0, v0, m, k, g, T, N, "euler");

t = linspace(0,T,N+1);

%PLOT 
figure;
ezplot(z,yAnalytic,[-0,T]); hold on;
plot(t,y); hold on;
legend('y (Analytic)','y (Numeric)'); title("Cromer | N=10000 | T=50");
xlabel('Time'); ylabel('Y-axis');

figure;
ezplot(z,vAnalytic,[-0,T]); hold on;
plot(t,v); hold on;
legend('v (Analytic)','v (Numeric)'); title("Cromer | N=10000 | T=50");
xlabel('Time'); ylabel('Velocity');

function [y,v] = NumericalFun(y0, v0, m, k, g, T, N, method)
    h = T/N;
    v = zeros(1,N+1);
    y = zeros(1,N+1);
    v(1) = v0;
    y(1) = y0;
    switch method
        case "euler";
            for j = 1:(N)
                v(j+1) = v(j) + h*(-k/m.*y(j));
                y(j+1) = y(j) + h*v(j);
            end
        case "cromer"
            for j = 1:(N)
                v(j+1) = v(j) + h*(-k/m.*y(j));
                y(j+1) = y(j) + h*v(j+1);
            end
        case "leapfrog"          
                vHalf = zeros(1,N);
                
            for j = 1:(N)
                vHalf(j) = v(j) + h*(-k/m*y(j))/2;
                y(j+1) = y(j) + h*vHalf(j);
                v(j+1) = vHalf(j) + h*((-k/m*y(j+1)))/2;
            end
        case "verlet"
            y(2) = y(1) + h*v(1);
            for j = 2:(N)
                y(j+1) = 2*y(j) - y(j-1) + h^2*(-k/m*y(j));
                v(j) = (y(j+1) - y(j-1))/(2*h);
            end
        otherwise
            fprintf("Try these for method:\n   euler\n   cromer\n   leapfrog\n   verlet\n");
    end
end
