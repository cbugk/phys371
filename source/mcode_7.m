%%%%% Assignment 7 - PHYS371 - 5th December 2018
%%%%% Karacan, Celil Bugra - 21401700

close all; clear; clc;

To = 100; % [K]
K = 237; % [W/mK]
C = 900; % [J/kgK]
rho = 2700; % [kg/m^3]
L = 1; % [m] - Length of wire
duration = 50; % [s] - Duration
h = 1; % there is no better guess, so why not one!

n = 31;
m = 3001;

dx = L/n;
dt = duration/m;

tempMatr = To*ones(m,n);
tempMatr(:,1) = 0;
tempMatr(:,n) = 0;
secndDiff = zeros(m,n);

tempMatrCF = tempMatr;
secndDiffCF = secndDiff;

tempMatrSum = tempMatr;
tempMatrSum2 = tempMatr;

tempMatrNwtn = tempMatr;
secndDiffNwtn = secndDiff;

stepSizeTitle = 'n = 31, m = 3001';

% PART-B-D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 2:m
    for j = 1:n-2
        secndDiff(i-1,j) = ( tempMatr(i-1,j+2) - 2*tempMatr(i-1,j+1) + tempMatr(i-1,j) )/dx^2;
        tempMatr(i,j) = tempMatr(i-1,j) + dt*K/C/rho*secndDiff(i-1,j);
    end
end

%plots - 1, 2
figure;
meshc(tempMatr);
xlabel('Location');
ylabel('Time');
title(stepSizeTitle);

figure;
contour(tempMatr);
title(stepSizeTitle);

% PART-E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 2:m
    for j = 2:n-1
        secndDiffCF(i-1,j) = ( tempMatrCF(i-1,j+1) - 2*tempMatrCF(i-1,j) + tempMatrCF(i-1,j-1) )/dx^2;
        tempMatrCF(i,j) = tempMatrCF(i-1,j) + dt*K/C/rho*secndDiffCF(i-1,j);
    end
end

%plots - 3, 4
figure;
meshc(tempMatrCF);
xlabel('Location');
ylabel('Time');
title(stepSizeTitle);

figure;
contour(tempMatrCF);
title(stepSizeTitle);

% PART-F
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% double loop for matrix evolves into triplke loop since there is a sum
% within.

for j = 1:n
    x  = dx*(j-1);
    for i = 1:m
        t = dt*(i-1);
        for b = 1:500
            tempMatrSum(i,j) = tempMatrSum(i,j) + (4*To/pi/b*sin(((2*b-1)*pi/L)*x)*exp(-((2*b-1)*pi/L)^2*K*t/C/rho));
        end
    end
end

%plots - 5, 6
figure;
meshc(tempMatrSum)
xlabel('Location');
ylabel('Time');
title(stepSizeTitle);

figure;
contour(tempMatrSum); 
title(stepSizeTitle);

% PART-G
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Combination of PART-E & PART-F


for i = 2:m
    x  = dx*(j-1);
    for j = 2:n-1
        t = dt*(i-1);
        for b = 1:500
            tempMatrSum2(i,j) = tempMatrSum2(i,j) + (4*To/pi/b*sin(((2*b-1)*pi/L)*x)*exp(-((2*b-1)*pi/L)^2*K*t/C/rho));
        end
        secndDiffNwtn(i-1,j) = ( tempMatrNwtn(i-1,j+1) - 2*tempMatrNwtn(i-1,j) + tempMatrNwtn(i-1,j-1) )*(n/L)^2;
        firstDer_t = K/C/rho*secndDiffNwtn(i-1,j) - h*tempMatrSum2(i,j);
        tempMatrNwtn(i,j) = tempMatrNwtn(i-1,j) + dt*firstDer_t;
    end
end

%plots - 7, 8
figure;
meshc(tempMatrNwtn);
xlabel('Location');
ylabel('Time');
title(stepSizeTitle);

figure;
contour(tempMatrNwtn);
title(stepSizeTitle);
