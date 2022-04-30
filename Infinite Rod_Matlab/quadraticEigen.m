function [K,X] = quadraticEigen(alpha1,alpha2,beta1,beta2,S0,N,R,n,rho,omega)

D11 = alpha2 + 2*alpha1;
D13 = alpha2 + beta2*S0;
D33 = alpha2 + 2*alpha1 + (1 + 4*beta1 + 2*beta2)*S0;
D55 = alpha1 + beta1*S0;
D66 = alpha1 + (1 + beta1)*S0;
dR = R/N;
% Annulus matrices and assemble
A = zeros(3*(N+1)); B = A; G = A; M = A;

for e = 1:N
    R1 = (e-1)*dR; R2 = e*dR;
    id = (1:6) + 3*(e-1)*ones(1,6);
    A(id,id) = A(id,id) + func_matrixA(D66,D33,R1,R2);
    B(id,id) = B(id,id) + func_matrixB(n,D13,D55,R1,R2);
    G(id,id) = G(id,id) + func_matrixG(n,alpha1,alpha2,D11,D55,R1,R2);
    M(id,id) = M(id,id) + func_matrixM(rho,R1,R2);
end

% Solve the quadratic eigenvalue problem
K = []; 
lap = size(omega, 2);
for i =1:lap
    A0 = G - omega(i)^2*M;
    A1 = 1i*B;
    A2 = A;
    [evec_quad,eval_quad] = polyeig(A0, A1, A2);
    % Choose wave numbers have negative imaginary parts
    indices = find((real(eval_quad)>0)&(imag(eval_quad)==0));
    K = [K; max(eval_quad(indices))];
    X = evec_quad(:, indices); % Eigenvector corresponding to omega(i)
end
