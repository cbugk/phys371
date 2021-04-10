%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Celil Bugra Karacan - 21401700 - PHYS371 Bilkent
%Homework-III Question-I Due: 29th Oct 2018
clear; clc;

syms T1 T2 T3 T4 L1 L2 L3 L4 M1 M2 M3 D theta1 theta2 theta3 theta4;

L1 = 3;    L2 = 2;     L3 = 4;     L4 = 5;  %[m]
M1 = 5;    M2 = 3;     M3 = 7;              %[kg]
D=8;                                        %[m]
g = 9.81;                                    %[m/s^2]

L = [L1 L2 L3 L4];
M = [M1; M2; M3];
%T = [T1; T2; T3; T4];
%%%%
A = [ 1 -1  0  0;
      0  1  1  0;
      0  0 -1  1];
  
A_pinv = pinv(A);

T_sin = A_pinv * M .* g;
%%%%
%K = [cos(theta1); cos(theta2); cos(theta3); cos(theta4)];
%D = L * K 
K = pinv(L).*D;
for i = 1:4;
    theta(i) = acos(K(i));
end

for j = 1:4
    test(j)=T_sin(j)*cot(theta(j));
end
test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Celil Bugra Karacan - 21401700 - PHYS371 Bilkent
%Homework-III Question-II Due: 29th Oct 2018
clear; clc;

syms R6; %for part two. Necessarry for ezplot().

%GIVEN CONSTANTS
R1=10; R2=20; R3=30; R4=40; R5=50; %R6=60;  %[Ohm]
V1=10; V2=20; V3=30; V4=40;                 %[V]

%MATRICES FORMED USING KIRCHOFF LAWS
M = [ R1  R2   0   0   0  R6;
       0   0   0  R4 -R5   0;
       0   0   0   0 -R5  R6;
       0   0   0  R4   0   0;
      R1   0  R3   0   0   0;
       1  -1  -1   0   0   0 ];

L = [  V1-V3-V4;
          0    ;
        V2-V3  ;
          -V4   ;
          V1   ;
          0     ];

% M * Current = L, then:     
Current = inv(M)*L

%P = I^2 * R 
Power6 = Current(6)^2 * R6

% Plotting Power6 symbolically
ezplot(Power6,[0,100]);
grid on; xlabel("R6 [ohms]"); ylabel("Power [Watt]");
