clear 
close all
clc
%%
syms  R R1 r D11 D12 D13 D33 D44 D55 D66 alpha1 alpha2 n rho 

N = [(R1-r)/(R1-R), 0, 0, (r-R)/(R1-R), 0, 0;
        0, (R1-r)/(R1-R), 0, 0, (r-R)/(R1-R), 0;
        0, 0, (R1-r)/(R1-R), 0, 0, (r-R)/(R1-R)];

dN = [-1/(R1-R), 0, 0, 1/(R1-R), 0, 0;
          0, -1/(R1-R), 0, 0, 1/(R1-R), 0;
          0, 0, -1/(R1-R), 0, 0, 1/(R1-R)];

% integral of A
A = N.'*diag([D66, D33, D66])*N*r;

% integral of B
B1 = dN.'*[0, -D13, 0; -D55, 0, 0; 0, 0, 0]*N*r; %B1

B2 = N.'*[0, D55, 0; D13, 0, 0; 0, 0, 0]*dN*r; %B2

B3 = N.'*[0, -D13, 0; D13, 0, -n*(D13+D55); 0, n*(D13+D55), 0]*N; %B3

% integral of G
G1 = dN.'*diag([D11, D55, alpha1])*dN*r; %G1

G2 = dN.'*[alpha2, 0, -n*alpha2; 0, 0, 0; n*alpha1, 0, -alpha1]*N; %G2

G3 = N.'*[alpha2, 0, n*alpha1; 0, 0, 0; -n*alpha2, 0, -alpha1]*dN; %G3

G4 = N.'*[n^2*alpha1+D11, 0, -n*(alpha1+D11); 0, n^2*D55, 0; -n*(alpha1+D11), 0, alpha1+n^2*D11]*N*1/r; %G4

% integral of M
M = N.'*diag([rho, rho, rho])*N*r;






