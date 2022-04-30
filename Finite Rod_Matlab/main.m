% Le Dinh Tan
% 12/06/2021
% Torsional spectrum for wave propagating in finite rods

clear, clc

%--------------------------------------------------------------------------
% Physical Parameters
%--------------------------------------------------------------------------
alpha1 = 1; % Shear modulus
beta1  = 0; % Material-dependent
S0 = 0*1e-3:2e-3:10e-3; % Initial stress 
rho = 1; % Density of steel
R = 0.0125; % Rod cross-section radius (m)
L = 1e3*R; % Rod length
N = 52; % The number of annuli 
dR = R/N;

for i=1:size(S0,2)
    %--------------------Annulus matrices and assemble---------------------
    A = zeros(N+1); G = A; M = A;
    for e = 1:N
        R1 = (e-1)*dR; R2 = e*dR;
        id = (1:2) + 1*(e-1)*ones(1,2); % id = (1:tong so bac tu do 2 node)+(so bac tu do cua 1 node)*(e-1)*ones(1,tong so bac tu do 2 node)
        A(id,id) = A(id,id) + func_matrixA(alpha1,beta1,S0(i),R1,R2);
        G(id,id) = G(id,id) + func_matrixG(alpha1,R1,R2,dR);
        M(id,id) = M(id,id) + func_matrixM(rho,R1,R2);
    end

    %------------Solve the quadratic eigenvalue problem--------------------
    omega = 0.001;
    % omega = 0.5;
    A0 = G - omega^2*M;
    A1 = zeros(N+1);
    A2 = A;
    [evec,eval] = polyeig(A0, A1, A2);

    %------Torsional spectrum (The torsional modes are antisymmetric)------
    indices = find((real(eval)>0)&(imag(eval)==0)|(imag(eval)<0));
    K = diag(eval(indices));
    X = evec(:, indices); 
    E = diag(exp(-1i*eval(indices)*L));

    %-----------------------Compute nodal forces---------------------------
    XE = [X, X*E; X*E, X];
    AX = [A*X*K, -A*X*K*E; -A*X*K*E, A*X*K];
    H = 2*pi*1i*AX/XE;
    U = [(0:N)'*dR; zeros(N+1,1)];
    P = H*U; % Nodal forces
    %-----------------------Compute the torque-----------------------------
    %***********************************************************
    % Load 1
    %***********************************************************
    K_theta = 0;
    for e=1:N+1
        K_theta = K_theta + P(e)*(e-1)*dR;
    end
    Iz = pi*R^4/2;
    K1(i) = K_theta/(alpha1*Iz/L);
    G_hat = alpha1*(1+(1+beta1)*S0(i)/alpha1);
    K2(i) = K_theta/(G_hat*Iz/L);
    %***********************************************************
    % Load 2
    %***********************************************************
    P1 = [zeros(N,1);1];
    H1 = H(1:N+1,1:N+1);
    U1 = H1\P1;
    K_theta1 = R^2/U1(N+1);
    K3(i) = K_theta1/(alpha1*Iz/L);
    K4(i) = K_theta1/(G_hat*Iz/L);
    K5(i) = K_theta1;
    K6(i) = K_theta;
end

for i=2:size(K5,2)
    K5(i) = K5(i)/K5(1);
end
K5(1) = 1;

for i=2:size(K5,2)
    K6(i) = K6(i)/K6(1);
end
K6(1) = 1;

Var1 = round(K1,3);
Var2 = round(K2,3);
Var3 = round(K6,3);
disp('Load 1')
table(S0', Var1', Var2', Var3','VariableNames', {'S0', 'Var1', 'Var2', 'Var3'})
disp('Load 2')
table(S0', round(K3,4)', round(K5,3)','VariableNames', {'S0', 'Var1', 'Var2'})




