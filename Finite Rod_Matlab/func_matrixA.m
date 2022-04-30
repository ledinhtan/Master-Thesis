function [A] = func_matrixA(alpha1,beta1,S0,R1,R2)

I1 = -R1^2/4 + (R1*R2)/6 + R2^2/12;
I2 = R1^2/12 - R2^2/12;
I3 = -R1^2/12 - (R1*R2)/6 + R2^2/4;
D1 = alpha1 + (1 + beta1)*S0;

A = [D1*I1 -D1*I2; -D1*I2 D1*I3];