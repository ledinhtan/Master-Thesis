% Le Dinh Tan
% 20/06/2021
% Torsional spectrum for wave propagating in infinite rods

clear 
close all
clc

% Physical Parameters
alpha1 = 1; % Shear modulus
alpha2 = 1; % Lame's constant
beta1  = -.5; % Material-dependent
beta2  = 0; % Initial-stress-dependent
S0 = .2*alpha1; % Initial stress
rho = .07; % Density 
R = 1;
n = 0; % Fourier number
N = 52; % The number of annuli 
Cs = sqrt(alpha1/rho); % The shear-wave velocity

% Neglect initial stress 
hold on
k = 0:.1:10;
Omega = sqrt((k*R).^2);
h = plot(k*R, Omega,'-b');
Omega1 = sqrt(5.136^2 + (k*R).^2);
plot(k*R, Omega1,'-b')
Omega2 = sqrt(8.417^2 + (k*R).^2);
plot(k*R, Omega2,'-b')


% The effect of initial stress 
omega = 1e-5:2*pi:10*Cs/R;
[K1,X] = quadraticEigen(alpha1,alpha2,beta1,beta2,S0,N,R,n,rho,omega);
k1 = abs(K1(:))*R;
Omega3 = sqrt(k1.^2);
h2 = plot(k1, Omega3,'--r');
Omega4 = sqrt(5.136^2 + k1.^2);
plot(k1, Omega4,'--r')
Omega5 = sqrt(8.417^2 + k1.^2);
plot(k1, Omega5,'--r')

xlabel('$kR$','Interpreter','latex','FontSize',14);
ylabel('$\frac{\omega R}{C_S}$','Interpreter','latex','FontSize',17)
title('Torsional Spectrum: beta1 = -1/2')
legend([h, h2],{'$S_0 = 0$', '$S_0 = 0.2G$'}, 'Interpreter', 'latex',...
            'FontSize', 17,...
            'location', 'SouthEast')
axis([0 10 0 10])


% print -deps nameoffigure % Do not have color when importing into latex
% print -depsc Figure5 % To have color when importing into latex

























