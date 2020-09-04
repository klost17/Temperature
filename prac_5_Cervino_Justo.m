%% Q1

% Critical value of K for which the Ra curve is minimum.
global Ra0
Ra0 = 1500;
nsKcr_even = fminbnd(@prac5criticalReEven,3,3.5);
nsRacr_even = prac5criticalReEven(nsKcr_even);

% Calculate the threshold heat. Coefficients corresponding to water.
g = 9.81; % Acceleration due to gravity in m s^-2.
d = 5e-3; % Separation between plates in m.
nu = 1.005e-6; % Kinematic viscosity in m^2 s^-1. At 20ºC.
kappa = 1.42e-7; % Thermal diffusivity in m^2 s^-1. At 20ºC.
alpha = 2.1e-4; % Thermal expansion in K^-1. At 20ºC.

Gamma = real(nsRacr_even*(kappa*nu)/(g*alpha*d^4));
DeltaTcr_even = Gamma*d

% For such small value of d, the threshold is of the order of 1.
% For higher values of d, the threshold is tiny, very close to 0.

%% Q2

Ra = real(nsRacr_even); K = real(nsKcr_even);
lambda = (Ra/K^4)^(1/3);
q1 = K*sqrt(1-lambda);
q2 = K*sqrt(lambda*(1+sqrt(3)*1i)/2+1);
q3 = K*sqrt(lambda*(1-sqrt(3)*1i)/2+1);

M(1,1)=cosh(q1/2);
M(2,1)=q1*sinh(q1/2);
M(3,1)=(q1^2-K^2)^2*(cosh(q1/2));
        
M(1,2)=cosh(q2/2);
M(2,2)=q2*sinh(q2/2);
M(3,2)=(q2^2-K^2)^2*(cosh(q2/2));
        
M(1,3)=cosh(q3/2);
M(2,3)=q3*sinh(q3/2);
M(3,3)=(q3^2-K^2)^2*(cosh(q3/2));

[V,~] = eig(M); % Columns of V are eigenvectors, and eigenvalues are in D.
A = V(1,2); % The second column of V corresponds to the eigenvector whose eigenvalue is 0.
B = V(2,2);
C = V(3,2);

% Define meshgrids for all plots, contourfs and quivers.
N = 51;
xmin = -1; xmax = 1;
[x,z] = meshgrid(linspace(xmin,xmax,N),linspace(-1/2,1/2,N));
[sf_x,sf_z] = meshgrid(linspace(xmin,xmax,N),linspace(-1/2,1/2,N));

% Compute
W_uppercase = A*cosh(q1*z) + B*cosh(q2*z) + C*cosh(q3*z);
w_hat = kappa/(Gamma*d^2)*W_uppercase;
D4 = q1^4*A*cosh(q1*z) + q2^4*B*cosh(q2*z) + q3^4*C*cosh(q3*z);
D2 = q1^2*A*cosh(q1*z) + q2^2*B*cosh(q2*z) + q3^2*C*cosh(q3*z);
T_hat = (1/(Ra*K^2))*(D4 - 2*K^2*D2 + K^4*W_uppercase);

% Normalization
normfact = max(max(abs(T_hat)));
T_hat = T_hat/normfact;
w_hat = w_hat/normfact;

Kn = K; Kx = K;
w_tilde = w_hat.*cos(Kn*x);
T_tilde = T_hat.*cos(Kx*x);

Dw_hat = kappa/(Gamma*d^2)*(q1*A*sinh(q1*z) + q2*B*sinh(q2*z) + q3*C*sinh(q3*z));
u_tilde = -(1/Kn)*Dw_hat.*sin(Kn*x);
u_tilde = u_tilde/normfact;

% Plotting
figure(1)
subplot(2,1,1)
plot(z(:,1),real(w_hat(:,1)));
hold on
subplot(2,1,2)
plot(z(:,1),real(T_hat(:,1)));
hold on

figure(2)
contourf(x,z,real(T_tilde)); colorbar;
hold on
quiver(x,z,real(u_tilde),real(w_tilde),'r')
hold off
xlabel('x'); ylabel('z'); title('1st mode (even solution) / No-slip BC'); grid

% Stress-free boundary conditions
n = 1;
sfKcr = n*pi/sqrt(2);
sfRacr = 27*(n*pi)^4/4;
Gamma = real(sfRacr*(kappa*nu)/(g*alpha*d^4));
sf_w_hat = kappa/(Gamma*d^2)*cos(n*pi*sf_z);
sf_T_hat = ((n*pi)^2+sfKcr^2)^2/(sfRacr*sfKcr^2)*cos(n*pi*sf_z);
normfact = max(max(sf_T_hat));
sf_T_hat = sf_T_hat/normfact;
sf_w_hat = sf_w_hat/normfact;
sf_w_tilde = sf_w_hat.*cos(sfKcr*sf_x);
sf_T_tilde = sf_T_hat.*cos(sfKcr*sf_x);
sf_u_tilde = n*pi*kappa/(Gamma*d^2*sfKcr)*sin(sfKcr*sf_x).*(sin(n*pi*sf_z));
sf_u_tilde = sf_u_tilde/normfact;

% Plotting
figure(1)
subplot(2,1,1)
plot(sf_z(:,1),sf_w_hat(:,1))
hold off
xlabel('z'); ylabel('w_{hat}(z)'); grid
title('1st mode (even solution) - Vertical velocity profile')
legend('No-slip BC', 'Stress-free BC')
subplot(2,1,2)
plot(sf_z(:,1),sf_T_hat(:,1))
hold off
xlabel('z'); ylabel('T_{hat}(z)'); grid
title('1st mode (even solution) - Temperature profile')
legend('No-slip BC', 'Stress-free BC')

figure(3)
contourf(sf_x,sf_z,sf_T_tilde); colorbar;
hold on
quiver(sf_x,sf_z,sf_u_tilde,sf_w_tilde,'r')
hold off
xlabel('x'); ylabel('z'); title('1st mode (even solution) / Stress-free BC'); grid

%% Q3a

% No-slip boundary conditions
Ra_ns_odd = 60000; ns_Kcr_odd = []; ns_Racr_odd = [];
for K = 2:0.01:8
    [solRa,it,res] = prac5newton(@prac5funOdd,Ra_ns_odd,K,1e-4,100);
    ns_Kcr_odd = [ns_Kcr_odd K];
    ns_Racr_odd= [ns_Racr_odd solRa];
end

% Critical value of K for which the Ra curve is minimum.
Ra0 = 20000;
nsKcr_odd = fminbnd(@prac5criticalReOdd,4,7);
nsRacr_odd = prac5criticalReOdd(nsKcr_odd);

% Plot everything
figure(4)
plot(ns_Kcr_odd,ns_Racr_odd,nsKcr_odd,nsRacr_odd,'ob');
xlabel('K'); ylabel('Ra_{cr}(K)'); grid
title('2nd mode (odd solution) - Neutral stability curve')
legend('No-slip BC')

%% Q3b & Q3c

Ra = real(nsRacr_odd); K = real(nsKcr_odd);
Gamma = real(nsRacr_odd*(kappa*nu)/(g*alpha*d^4));
lambda = (Ra/K^4)^(1/3);
q1 = K*sqrt(1-lambda);
q2 = K*sqrt(lambda*(1+sqrt(3)*1i)/2+1);
q3 = K*sqrt(lambda*(1-sqrt(3)*1i)/2+1);

M(1,1)=sinh(q1/2);
M(2,1)=q1*cosh(q1/2);
M(3,1)=(q1^2-K^2)^2*(sinh(q1/2));
        
M(1,2)=sinh(q2/2);
M(2,2)=q2*cosh(q2/2);
M(3,2)=(q2^2-K^2)^2*(sinh(q2/2));
        
M(1,3)=sinh(q3/2);
M(2,3)=q3*cosh(q3/2);
M(3,3)=(q3^2-K^2)^2*(sinh(q3/2));

[V,D] = eig(M); % Columns of V are eigenvectors, and eigenvalues are in D.
A = V(1,2); % The second column of V corresponds to the eigenvector whose eigenvalue is 0.
B = V(2,2);
C = V(3,2);

% Compute
W_uppercase = A*sinh(q1*z) + B*sinh(q2*z) + C*sinh(q3*z);
w_hat = kappa/(Gamma*d^2)*W_uppercase;
D4 = q1^4*A*sinh(q1*z) + q2^4*B*sinh(q2*z) + q3^4*C*sinh(q3*z);
D2 = q1^2*A*sinh(q1*z) + q2^2*B*sinh(q2*z) + q3^2*C*sinh(q3*z);
T_hat = (1/(Ra*K^2))*(D4 - 2*K^2*D2 + K^4*W_uppercase);
w_hat = imag(w_hat);  % NO SÉ POR QUÉ PERO SON IMAGINARIOS PUROS
T_hat = imag(T_hat);  % NO SÉ POR QUÉ PERO SON IMAGINARIOS PUROS

% Normalization
normfact = max(max(abs(T_hat)));
T_hat = T_hat/normfact;
w_hat = w_hat/normfact;

Kn = K; Kx = K;
w_tilde = w_hat.*cos(Kn*x);
T_tilde = T_hat.*cos(Kx*x);

Dw_hat = kappa/(Gamma*d^2)*(q1*A*cosh(q1*z) + q2*B*cosh(q2*z) + q3*C*cosh(q3*z));
u_tilde = -(1/Kn)*Dw_hat.*sin(Kn*x);
u_tilde = imag(u_tilde);   % NO SÉ POR QUÉ PERO SON IMAGINARIOS PUROS
u_tilde = u_tilde/normfact;

% Plotting
figure(5)
subplot(2,1,1)
plot(z(:,1),real(w_hat(:,1)))
hold on
subplot(2,1,2)
plot(z(:,1),real(T_hat(:,1)))
hold on

figure(6)
contourf(x,z,real(T_tilde)); colorbar;
hold on
quiver(x,z,real(u_tilde),real(w_tilde),'r')
hold off
xlabel('x'); ylabel('z'); title('2nd mode (odd solution) / No-slip BC'); grid

% Stress-free boundary conditions
n = 2;
sfKcr = n*pi/sqrt(2);
sfRacr = 27*(n*pi)^4/4;
Gamma = real(sfRacr*(kappa*nu)/(g*alpha*d^4));
sf_w_hat = kappa/(Gamma*d^2)*sin(n*pi*sf_z);
sf_T_hat = ((n*pi)^2+sfKcr^2)^2/(sfRacr*sfKcr^2)*sin(n*pi*sf_z);
normfact = max(max(sf_T_hat));
sf_T_hat = sf_T_hat/normfact;
sf_w_hat = sf_w_hat/normfact;
sf_w_tilde = sf_w_hat.*cos(sfKcr*sf_x);
sf_T_tilde = sf_T_hat.*cos(sfKcr*sf_x);
sf_u_tilde = -n*pi*kappa/(Gamma*d^2*sfKcr)*sin(sfKcr*sf_x).*(cos(n*pi*sf_z));
sf_u_tilde = sf_u_tilde/normfact;

% Plotting
figure(5)
subplot(2,1,1)
plot(sf_z(:,1),sf_w_hat(:,1))
hold off
xlabel('z'); ylabel('w_{hat}(z)'); grid
title('2nd mode (odd solution) - Vertical velocity profile')
legend('No-slip BC', 'Stress-free BC')
subplot(2,1,2)
plot(sf_z(:,1),sf_T_hat(:,1))
hold off
xlabel('z'); ylabel('T_{hat}(z)'); grid
title('2nd mode (odd solution) - Temperature profile')
legend('No-slip BC', 'Stress-free BC')

figure(7)
contourf(sf_x,sf_z,sf_T_tilde); colorbar;
hold on
quiver(sf_x,sf_z,sf_u_tilde,sf_w_tilde,'r')
hold off
xlabel('x'); ylabel('z'); title('2nd mode (odd solution) / Stress-free BC'); grid

%% Q4

% Analytical