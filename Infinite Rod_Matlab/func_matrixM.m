function [M] = func_matrixM(rho,R1,R2)

I1 = -R1^2/4 + (R1*R2)/6 + R2^2/12;
I2 = 1/12*(R1-R2)*(R1+R2);
I3 = -R1^2/12 - (R1*R2)/6 + R2^2/4;
Ro = eye(3)*rho;

M  = [I1*Ro -I2*Ro; -I2*Ro I3*Ro];

