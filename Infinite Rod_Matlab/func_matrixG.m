function [G] = func_matrixG(n,alpha1,alpha2,D11,D55,R1,R2)

% Integrals
I1 = (R2 + R1)/(2*(R2 - R1));
I2 = 1/2;
I3 = -1/2;

I4 = @(r) [((alpha1*n^2 + D11)*(R2 - r)^2)/(r*(R2 - R1)^2), 0, -(n*(D11 + alpha1)*(R2 - r)^2)/(r*(R2 - R1)^2),...
          -((alpha1*n^2 + D11)*(R1 - r)*(R2 - r))/(r*(R2 - R1)^2), 0, (n*(D11 + alpha1)*(R1 - r)*(R2 - r))/(r*(R2 - R1)^2);...
          0, (D55*n^2*(R2 - r)^2)/(r*(R2 - R1)^2), 0, 0, -(D55*n^2*(R1 - r)*(R2 - r))/(r*(R1 - R2)^2), 0;...
          -(n*(D11 + alpha1)*(R2 - r)^2)/(r*(R1 - R2)^2), 0, ((D11*n^2 + alpha1)*(R2 - r)^2)/(r*(R1 - R2)^2),...
          (n*(D11 + alpha1)*(R1 - r)*(R2 - r))/(r*(R1 - R2)^2), 0, -((D11*n^2 + alpha1)*(R1 - r)*(R2 - r))/(r*(R1 - R2)^2);...
          -((alpha1*n^2 + D11)*(R1 - r)*(R2 - r))/(r*(R1 - R2)^2), 0, (n*(D11 + alpha1)*(R1 - r)*(R2 - r))/(r*(R1 - R2)^2),...
          ((alpha1*n^2 + D11)*(R1 - r)^2)/(r*(R1 - R2)^2), 0, -(n*(D11 + alpha1)*(R1 - r)^2)/(r*(R1 - R2)^2);...
          0, -(D55*n^2*(R1 - r)*(R2 - r))/(r*(R1 - R2)^2), 0, 0, (D55*n^2*(R1 - r)^2)/(r*(R1 - R2)^2), 0;...
          (n*(D11 + alpha1)*(R1 - r)*(R2 - r))/(r*(R1 - R2)^2), 0, -((D11*n^2 + alpha1)*(R1 - r)*(R2 - r))/(r*(R1 - R2)^2),...
          -(n*(D11 + alpha1)*(R1 - r)^2)/(r*(R1 - R2)^2), 0, ((D11*n^2 + alpha1)*(R1 - r)^2)/(r*(R1 - R2)^2)]; 
      

% D matrices
Dg1 = diag([D11, D55, alpha1]);
Dg2 = diag([alpha2, 0, -alpha1]); 
Dg2(1,3) = -n*alpha2;
Dg2(3,1) =  n*alpha1;
Dg3 = diag([n^2*alpha1+D11, n^2*D55, alpha1+n^2*D11]);
Dg3(1,3) = -n*(alpha1+D11);
Dg3(3,1) = Dg3(1,3);

% G components and G
G1 = I1*[Dg1    -Dg1;   -Dg1    Dg1];
G2 = [-I2*Dg2 I3*Dg2; I2*Dg2 -I3*Dg2];
G3 = G2';
if R1 == 0
    G4 = I4(R2/2)*R2;
else
    I5 =  R2^2*(log(R2)-log(R1))/(R1-R2)^2 - (R1-3*R2)/(2*(R1-R2));
    I6 = (R1+R2)/(2*(R1-R2)) + (2*R1*R2*(log(R2)-log(R1)))/(2*(R1-R2)^2);
    I7 = (3*R1-R2)/(2*(R1-R2)) + (R1^2*(log(R2)-log(R1)))/(R1-R2)^2;
    G4 = [Dg3*I5, Dg3*(-I6); Dg3*(-I6), Dg3*I7];
end

G  = G1 + G2 + G3 + G4;

