close all; clear; clc;

M = [25 50 75 100 125 150 175 200 225 250];
n = 50;

for h = 1:length(M)
m = M(h)
    for j = 1:n % 'n' many trials per a given 'm'    
        t = 0; % iterations passed
        w1 = ceil(m*rand());
        w2 = w1;
        while w1 == w2
            w2 = ceil(m*rand()); % w1 ~= w2, before evaluation
        end
        Distance = w2 - w1; % initial distance (t = 0)
        walker = [ w1, w2 ]; % walkers spawned
                             % (easier to use with 'for loop's)
        Journey = walker;  % Positions of walkers are logged
        while walker(1) ~= walker(2)
            t = t+1;
            b = ceil(2*rand()); % creating a half probability |*|
            for i = 1:2  % per walker                         |*|               
                % if statements checking for boundaries
                switch b
                    case 1 
                        if walker(i) ~= 1             
                            walker(i) = walker(i) -1;
                        end 
                    case 2
                        if walker(i) ~= m
                            walker(i) = walker(i) +1;
                        end
                end
            end    
            Journey(t+1, 1) = walker(1);
            Journey(t+1, 2) = walker(2);
            d = abs(walker(1) - walker(2));
            Distance = [Distance d];
        end
        % max(Distance) rather than abs(w1 - w2),
        % so the code is applicable to 'drunk' walkers.
        maxD = max(Distance);
        T(j) = t; 
        time = 0:t;
    
        %figure;
        %subplot(1,2,1);
        %plot(Journey(:,1),time,':r',Journey(:,2),time,':b')
        %xlabel('Position');
        %ylabel('Iteration');
        %axis([0, m+1, 0, t])
        %subplot(1,2,2)
        %scatter(Distance,time,'.')
        %xlabel('Distance');
        %ylabel('Iteration');
        %axis([0, maxD, 0, t])
        %suptitle(['m = ' num2str(m) ', t = ' num2str(t)])
    end
    
    meanT(h) = sum(T)/n;
    for p = 1:n
        sqDiff(p) = (T(p)-meanT(h))^2;
    end
    err(h) = sqrt(sum(sqDiff) /n);
    
end

figure;
plot(err)
xlabel('Members of M (length set)');
ylabel('Standard Deviation');
title(['n = ' num2str(n)])

figure;
errorbar(M,meanT,err,'r');
xlabel('Length of Finite Line');
ylabel('Mean Value of Iterations');
title(['n = ' num2str(n)])
%errorbarlogx;
%fucntion by Frederic Moisy (above) can be obtained online.
%in case of a logarithmic increase in 'm'.