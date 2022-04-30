function [B] = func_matrixB(n,D13,D55,R1,R2)

% Integrals
a1 = R1 + R2;
a2 = R2 - R1;
I4 = (R2*a1 - 2*R1^2)/(6*a2);
I5 = (R1*a1 - 2*R2^2)/(6*a2);
I6 = a2/3;
I7 =  -a2/6;

% D matrices
Db1 = zeros(3); Db1(1,2) = -D13; Db1(2,1) = -D55;
Db2 = -Db1';
Db3 = zeros(3); 
Db3(1,2) = -D13; Db3(2,1) = D13;
Db3(2,3) = -n*(D13 + D55);
Db3(3,2) =  n*(D13 + D55);

% B components and B
B1 = [-I4*Db1 I5*Db1; I4*Db1 -I5*Db1];
B2 = [-I4*Db2 I4*Db2; I5*Db2 -I5*Db2];
B3 = [I6*Db3  -I7*Db3; -I7*Db3 I6*Db3];

B  = B1 + B2 + B3;

