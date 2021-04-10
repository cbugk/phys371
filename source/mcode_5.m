%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Celil Bugra Karacan - 21401700 - PHYS371 Bilkent
%Homework-V Question-I Due: 12th Nov 2018
clear; close all; clc;


h = 1/100; % step size.
m=0; R=0;

% Non-variable Coefs.
rho = 1.2; % density of atmosphere (25 deg. Celcius). [kg/meter^3]
g = 9.80; % gravitational acceleration [meter/s^2].

global N;% globalized to be implemented through functions.
balls = ["basketball", "baseball", "bowlingball"];

b=3; %% when for loop below is commented out sets ball type.
%for b = 1:3 % (**^^)

    % Variable Coefs.
    [m R] = setBall(balls(b)); % sets "A" and "m" variables globally.
    Cd = 0.5; % drag coefficient. (given as D in question).
    A = pi*R^2; % cross section of sphere [meter^2]
    
    FGC = rho*Cd*A/2; % F_air Global Coef, in f_air func.
    ARC = 2*m/(Cd*rho*A); % Analytical Repeated Coef, for readablity. 
    
    % initial conditions
    y0 = 0;
    v0 = 0;
    
    
    
    %%% NUMERICAL SOLN. {
    [yNum vNum] = NumericalFun(h, m, g, y0, v0,FGC);
    N = length(vNum); % Number of steps
    T = h * N; % Total time enlapsed [seconds]
    fprintf("Last Position: %f\nNumber of steps: %f\nTime enlapsed: %f\n",yNum(length(yNum)),N,T);
    %%% }
    
    %%% ANALYTICAL SOLN. {
    [yAn vAn] = AnalyticalFun(ARC,N,h,g);
    %%% }
    
    
    %%% PLOTTING {
    t = [0];  % linspace() required integers, so for loop.
    %for i = 1:N-1 (use instead so ball stops at y=440) (***)
    for i = 1:N-1
    %for i = 1:1591 %(h = 1/100)
        t(i+1) = i*h;
    end
    
    figure(1); title("No Air Resistance");
    plot(t,yAn); hold on
    plot(t,yNum); hold on
    legend('Location','NorthWest'); %Comment out legend() down below and
                            %remove comment from outermost for() loop  (**^^)
    %legend("y(Analytical)","y(Numerical)",'Location','NorthWest'); 
    xlabel('Time'); ylabel('Y-axis'); grid minor;
    
    figure(2); title("No Air Resistance");
    plot(t,vAn); hold on;
    plot(t,vNum); hold on;
    legend('Location','NorthWest'); %Comment out legend() down below and
                           %remove comment from outermost for() loop  (**^^)
    %legend("v(Analytical)","v(Numerical)",'Location','NorthWest');
    xlabel('Time'); ylabel('Velocity'); grid minor;
    %%% }
%end   % (**^^)

    
function [m R] = setBall(ball)
    m=0; R=0;% respectively mass & radius
    switch ball
        case "basketball"
            R = 0.012;
            m = 0.620;
        case "baseball"
            R = 0.038;
            m = 0.145;
        case "bowlingball"
            R = 0.01083;
            m = 7.200;
        otherwise
            fprintf("Try these balls:\n   basketball\n   baseball\n   bowlingball\n");
    end
end

function [f] = f_air (vel, FGC)
    f = FGC *vel^2;
end
 
function [y,v] = NumericalFun(h, m, g, y0, v0,FGC)
    v = [v0];
    y = [y0];
    
    j=1; % will end up as size(v)
    while   y(j) <= 440 %(use instead so ball stops at y=440)  (***)
    %while  j <= 1591 %(**^^)
        v(j+1) = v(j) + h*(g - f_air(v(j),FGC)/m);
        y(j+1) = y(j) + h*v(j);
        j = j+1;
    end
end

function [y,v] = AnalyticalFun(ARC,N,h,g)
    y = [];
    v = [];
    for j = 1:N
        y(j) = ARC * log(cosh((h*(j-1)) * sqrt(g/ARC)));
        v(j) = sqrt(ARC*g) * tanh( sqrt((h*(j-1)) * sqrt(g/ARC)));
    end
end