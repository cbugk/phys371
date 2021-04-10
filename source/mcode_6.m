%%%% Karacan Celil Bugra 21401700 26Nov18 Phys371

clear; close all; clc;
%%%%

n=250;
t = linspace(0,2*pi,n); % range of period
s = 1./(1 - 0.9*sin(t)); % func given s(t)
alpha = 5; % for alpha = 0, y = s
y = s;
for j = 1:n
    y(j) = s(j) + alpha*rand(); %noice added
end
%%%%


S = fft(y,n); % Fourier transform S(omega)

S_Abs_Sq = abs(S).^2;


%{
%% CODE BELLOW takes too long to run, and performs same (for n=3) as autocorr()
%% so it's not used yet will remain as a proof of work.

syms k;
tau = t;
z1 = 1./(1 - 0.9*sin(k));
acf = zeros(1,n);
tic;                    %to measure how long it takes
for i = 1:n
    acf(i) = int(z1*(1./(1 - 0.9*sin(k + tau(i)))), k, 0, 2*pi);
end
toc;
%%%
%also xcorr() can be used instead of autocorr() but
%range2 = linspace(0,2*pi,2*n-1); % Since xcorr() gives 2n-1 terms
%}

% ,NumLags, n-1) is necessary to make sizes of t and acf the same
acf = autocorr(s,'NumLags',n-1); %part-c
A = fft(acf);

%%%%

myTitle = '\alpha = 5'

figure;
plot(t,y);
xlim([-0.05 0.05+2*pi]);
legend('s(t)');
xlabel('Time [s]');
title(myTitle);


figure; % part-a
plot(t,S);
xlim([-0.05 0.05+2*pi]);
legend('S(\omega)');
xlabel('Time [s]');
title(myTitle);

%{
figure; %part-b
semilogy(t,S_Abs_Sq);
xlim([-0.05 0.05+2*pi]);
legend('|S(\omega)|^2');
xlabel('Time [s]');
%}

figure; % part-b and part-d
semilogy(t,S_Abs_Sq); hold on;
semilogy(t,A/sqrt(2*pi));
legend('|S(\omega)|^2','A(\omega) /sqrt(2\pi)')
xlim([-0.05 0.05+2*pi]);
xlabel('Time [s]');
title(myTitle);

figure; % part-d
plot(t,A);
legend('A(\omega)')
xlim([-0.05 0.05+2*pi]);
xlabel('Time [s]');
title(myTitle);



