clc
clear
close all

a1=sqrt(2);b1=1;K=-1;
C1=10^-9; fc=2*10^3;

% C2=C1*4*b1*(1-K)/(a1^2);
C2=10^-8;
R2= (a1*C2 - ((a1^2) * (C2^2) - (4*b1*C1*C2*(1-K)))^0.5)/(4*pi*fc*C1*C2);
R1 = -R2/K;
R3= b1/(4*(pi^2)*(fc^2)*C1*C2*R2);
