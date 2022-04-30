function [A] = func_matrixA(D66,D33,R1,R2)

I1 = -R1^2/4 + (R1*R2)/6 + R2^2/12;
I2 = -1/12*(R1 - R2)*(R1 + R2);
I3 = -R1^2/12 - (R1*R2)/6 + R2^2/4;
Da = diag([D66,D33,D66]);

A  = [I1*Da I2*Da; I2*Da I3*Da];

