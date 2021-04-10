%Q1
% CELIL BUGRA KARACAN - PHYS371 - ASSIGN-2
clear; clc;

%Given constants
a0 = 5.2917721067*10^(-11); % Bohr Radius [meter]
a = 20 * a0;
v0 = 13.6; % [eV]
hBar = 1.05457148*10^(-34); % [meter^2 *kg /s]
eMass = 9.10938356*10^(-31); % Electron Mass [kg]

% Equations given
% eqnEven = sqrt(E + v0) .* tan( a./hBar.*sqrt(2.*eMass.*(E + v0) ) ) - sqrt(-E);
% eqnOdd  = -sqrt(E + v0) .* cot( a./hBar.*sqrt(2.*eMass.*(E + v0) ) ) - sqrt(-E);

% Array of eigen values initiated.
E_Array = [-14 -14 -14 -14 -14 -14 -14 -14 -14 -14];

%Linspace of allowed energy levels with as much points as possible (w/out jepordizing run time) 
E = linspace(-13.6,0,136000001);



c=0;

for i = 1:length(E)    
    eqnEven = sqrt(E(i) + v0) .* tan( a./hBar.*sqrt(2.*eMass.*(E(i) + v0) ) ) - sqrt(-E(i));
    if -0.0000009 <= eqnEven  && eqnEven <= 0.0000009
       c = c + 1;
       if c <= length(E_Array)
            for j = 1:length(E_Array)
                if E(i) > E_Array(j)
                    E_Array(c) = E(i);
                    break
                end
            end
       end
    end       
end

fprintf("-----------------------\n");
E_Array
fprintf("-----------------------\n");


z = length(E_Array) - c;

m=0;
n=0;

for i = 1:length(E)
    
    eqnOdd = sqrt(E(i) + v0) .* tan( a./hBar.*sqrt(2.*eMass.*(E(i) + v0) ) ) - sqrt(-E(i));
    if -0.0000009 <= eqnOdd  && eqnOdd <= 0.0000009
       n = n + 1;
       if n <= length(E_Array)
            for j = n+1:c
                if E(i) > E_Array(j) && m < z
                    m = m +1;
                    E_Array(c+m) = E(i);
                    break
                end
            end
       end
    end       
end

fprintf('%f energy levels are determined\n',c+m);
E_Array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Q2
% CELIL BUGRA KARACAN - PHYS371 - ASSIGN-2
clear; clc;

%defining symvars x and y
syms x y;

%Declaring constant given
L1 = 0.1;  % [m]
L2 = 0.1;  % [m]
D  = 0.1;  % [m]
m  = 0.1;  % [kg]
g  = 9.81; % [m/s^2]
k1 = 10;   % [N/m]

%Part-A set k2
k2 = 20;   % [N/m]

%Part-B creating vector k2, and same-length x, y, V values.
%k2 = linspace(0,20,51);  %%uncomment for Part-B
solnX = zeros(1,length(k2));
solnY = solnX;
V = solnX;


% for loop to execute each k2 value and find correspanding x, y, V values.
for i = 1:length(k2)       
    
    %given fuction for Potential
    V_fun = @(x,y) 0.5 .* k1 .* (sqrt(x^2+y^2) - L1)^2 + 0.5 .* k2(i) .*(sqrt((x - D)^2 + y^2) - L2)^2 - m.*g.*y;
    
    %taking gradient of V_fun with respect to (x,y)
    gra = -1 * gradient(V_fun(x,y));
    
    %solving for gra == 0 condition by vpasolve. ( Matlab suggested after trying solve() )
    soln = vpasolve(gra,[x y]);
    
    solnX(i) = soln.x; % satisfying x value
    solnY(i) = soln.y; % satisfying y value
end

fprintf('Part-A\n(x,y) = (%f, %f)',solnX(length(k2)),solnY(length(k2)));


%plot(solnX,solnY,'.')
%grid on; xlabel('x [m]'); ylabel('y [m]');
%hold on; legend('k2 C [0,20] N/m', 'k2 = 20 N/m');

%Part-B
plot(k2, solnX, '+')
hold on;
plot(k2, solnY, '.')
grid on; xlabel('k2 [N/m]');
legend('x', 'y');
