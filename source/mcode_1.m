%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%karacan_2018_Oct_15_hw1_q1.m
%Celil Bugra Karacan - 21401700 - PHYS371 Bilkent
%Homework-I Question-1  Due: 15th Oct 2018

clear; clc;
syms theta; %defining theta as main symbolic variable.
eps=[0 0.1 0.5 0.9]; %vector of values desired for epsilon.
labels=["eps=0.0","eps=0.1","eps=0.5","eps=0.9"]; %for legend of plots.
TitL = ["(a)","","",""]; TitR = ["(b)","","",""]; %not recommended way of\
                                                  %adding two column titles
for i=1:4; %running for each value of eps
    r(theta) = (1-eps(i).^2)/(1-eps(i).*cos(theta)); %given by the question
    x = r.*cos(theta); %converting from polar to cartesian coordinates
    y = r.*sin(theta); %both symbolical functions of theta as r.

    subplot(4,2,2*i-1); %plots on left column
    ezplot(r,[0 2*pi]); %{r versus theta} plotting symbolically
    title(TitL(i)); %title visible at top left
    legend(labels(i)); grid on;       
    xlabel('theta'); ylabel('r(theta)');

    subplot(4,2,2*i);  %plots on right column
    ezplot(x,y); %{x versus y} plotting symbolically
    title(TitR(i)); %title visible at top right
    legend(labels(i)); grid on;    
    xlabel('x = r*cos(theta)'); ylabel('y = r*sin(theta)');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%karacan_2018_Oct_15_hw1_q2.m
%Celil Bugra Karacan - 21401700 - PHYS371 Bilkent
%Homework-I Due: 15th Oct 2018
N = 60;
x = [10 2 -2 -10];
labels = ["x = 10", "x = 2", "x = -2", "x = -10"]; %legend labels for plots
TitL = ["err\_plot", "", "", ""];  %not recommended way of adding 
TitR = ["err\_plot\_mod", "", "", ""];         %two column titles

for j=1:length(x)   %for loop to plot all values in a figure
    subplot(4,2,2*j-1);
    err_plot(x(j),N); grid on;
    title(TitL(j)); %only TitL(1) is visible (at top left)
    legend(labels(j)); xlabel('N'); ylabel('absolute fractional error');
    
    subplot(4,2,2*j)
    err_plot_mod(x(j),N); grid on;
    title(TitR(j)); %only TitR(1) is visible (at top right)
    legend(labels(j)); xlabel('N'); ylabel('absolute fractional error');
end

function [] = err_plot(x,N) %Part-A 
    err = zeros(1,N+1); %length of N+1
    n = linspace(0,N,N+1); %length of N+1
    sum = 1; %S(x,N) = 1, independent from x.
    for i=0:N
        sum = sum + (x.^i)./factorial(i); %incrementing partial sum 
        err(i+1) = abs(sum - exp(x)) ./ exp(x); %storing corresponding err
    end
    plot(n,err,'r'); %void function, only plots N vs err.
end

function [] = err_plot_mod(x,N) %Part-B
    err = zeros(1,N+1);
    n = linspace(0,N,N+1);
    sum = 1;
    if x>=0
        for i=0:N
            sum = sum + (x.^i)./factorial(i);
            err(i+1) = abs(sum - exp(x)) ./ exp(x);
        end
    else %for x<0
        for i=0:N
            sum = sum + (-x).^i./factorial(i);
            err(i+1) = abs(1./sum - exp(x)) ./ exp(x); %different approx     
        end                                            %is used for e^x :
    end                                     % 1/S(-x,N) = 1/e^(-x) = e^x
    plot(n,err,'b'); %again void function that plots.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%