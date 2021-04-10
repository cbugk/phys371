%%% CELIL BUGRA KARACAN * 21401700 * 2018-XII-20 assign9
%%% PHYS 371

close all; clear; clc;


k = 15; % number of trials per a given 's'
% Array of different s'es are below
throwsArr = [ 10^3, 5*10^3, 10^4, 5*10^4, 10^5, 5*10^5, 10^6, 5*10^6, 10^7, 5*10^7 ]

for m = 1:length(throwsArr)
    s = throwsArr(m); % number of total shots (all fallin the square)
    numPi = []; % array of approximations for pi (length(numPi)=k) 

    for j = 1:k;
        c = 0; % initializing number of hits within the circle
        for i = 1:s
            x = rand();
            y = rand();
            if x^2 + y^2 < 1
                c = c +1; % incrementing
            end
        end
        numPi(j) = 4 *c /s;
    end
    meanPi(m) = sum(numPi) /k;
    
    for n = 1:k
        variance(n) = (numPi(n)-meanPi(m))^2;
    end
    err(m) = sqrt(sum(variance) /k);
end
figure;
plot(err)
xlabel('Trials');
ylabel('Standard Deviation');

figure;
errorbar(throwsArr,meanPi,err,'r');
xlabel('# of Throws');
ylabel('Mean Value For \pi');
title('errorbar()');

figure;
errorbar(throwsArr,meanPi,err,'r');
xlabel('# of Throws');
ylabel('Mean Value For \pi');
title('errorbarlogx()');
errorbarlogx;







%%% Below is Frederic Moisy's code for plotting errorbar with semilogx kind
%%% of manner. code below is taken from address below, for esthetics
%%% reasons
%%% https://www.mathworks.com/matlabcentral/fileexchange/9715-errorbarlogx-m

    function errorbarlogx(epsilon)
    %ERRORBARLOGX Homogenize the error bars for X-axis in log scale.
    %   ERRORBARLOGX turns the X-axis of the current error bar plot to log
    %   scale, and homogonizes the length of the horizontal segements which
    %   terminate the vertical error bars.
    %
    %   By default, Matlab's ERRORBAR draws vertical error bars which are
    %   terminated by small horizontal segments of uniform length for the X-
    %   axis in linear scale. But when turning the X-axis to log scale, these
    %   segments become uneven. Using ERRORBARLOGX makes them uniform again.
    %
    %   ERRORBARLOGX(N) specifies the relative length of the horizontal
    %   segments, normalized with the total range of the data. By default,
    %   N=0.01 is used.
    %
    %   Limitations: ERRORBARLOGX acts only on the last drawn curve. If this
    %   curve is not an error bar plot, it won't work.
    %
    %   Example:
    %      x=logspace(1,3,20);
    %      y=5*(1 + 0.5*(rand(1,20)-0.5)).*x.^(-2);
    %      errorbar(x,y,y/2,'o-');
    %      errorbarlogx(0.03);
    %       
    %   F. Moisy
    %   Revision: 1.00,  Date: 2006/01/20
    %
    %   See also ERRORBAR.
    % History:
    % 2006/01/20: v1.00, first version.
    if nargin==0, epsilon=0.01; end; % default normalized segment length
    
    set(gca,'XScale','log'); % set the X axis in log scale.
    ca = get(gca); % current axe properties
    cd = get(ca.Children(1)); % current data properties
    heb = cd.Children(2); % handle to current error bars
    ceb = get(heb); % current error bars properties
    % ceb.XData is an array of length 9*length(data).
    % for each data point, ceb.XData contains 3 blocks of 3 numbers:
    % X0 X0 NaN X0-DX X0+DX Nan XO-DX X0+DX Nan
    % So it is necessary to change the 4th, 5th, 7th and 8th values.
    % logarithmic length of the horizontal segments:
    dx = 10^(log10(ceb.XData(length(ceb.XData)-8)/ceb.XData(1))*epsilon);
    % computes the new horizontal segments for each data point:
    for i=1:(length(ceb.XData)/9),
        ii=(i-1)*9+1; % index of the first error bar data for the point #i
        ceb.XData(ii+3) = ceb.XData(ii)/dx;
        ceb.XData(ii+4) = ceb.XData(ii)*dx;
        ceb.XData(ii+6) = ceb.XData(ii)/dx;
        ceb.XData(ii+7) = ceb.XData(ii)*dx;
    end;
    set(heb,'XData',ceb.XData); % sets the changes in the error bar properties
    end
    
 %%%
