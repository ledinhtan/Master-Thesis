function [M] = func_matrixM(rho,R1,R2)

I1 = -R1^2/4 + (R1*R2)/6 + R2^2/12;
I2 = -R1^2/12 + R2^2/12;
I3 = -R1^2/12 - (R1*R2)/6 + R2^2/4;

M  = [I1*rho I2*rho; I2*rho I3*rho];
