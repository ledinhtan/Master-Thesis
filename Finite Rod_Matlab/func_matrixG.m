function [G] = func_matrixG(alpha1,R1,R2,dR)

% Integrals
I1 = (R2^2/2 - R1^2/2)/((R1 - R2)^2);
I2 = R1^2/(2*(R1 - R2)^2) - R2/(R1 - R2) - R2^2/(2*(R1 - R2)^2);
I3 = R1^2/(2*(R1 - R2)^2) - R1/(R1 - R2) - R2^2/(2*(R1 - R2)^2);

I4 = @(r) [(alpha1*(R2 - r)^2)/(r*(R1 - R2)^2), -(alpha1*(R1 - r)*(R2 - r))/(r*(R1 - R2)^2);
          -(alpha1*(R1 - r)*(R2 - r))/(r*(R1 - R2)^2),(alpha1*(R1 - r)^2)/(r*(R1 - R2)^2)];

% G components and G
G1 = [alpha1*I1    -alpha1*I1;   -alpha1*I1    alpha1*I1];
G2 = [alpha1*I2 -alpha1*I3; -alpha1*I2 alpha1*I3];
G3 = G2';

%=======================================================
% Midpoint rule to neglect singularity
%=======================================================
N = 20;% The number of subinterval each annuli
dr = dR/N;
G4 = zeros(2);
for i=1:N
    midpoint = R1 + i*dr - dr/2;
    G4 = G4 + I4(midpoint)*dr;
end
G  = G1 + G2 + G3 + G4;

