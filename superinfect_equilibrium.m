% equilibrium for superinfection parameter setting

D=0.005; r=0.1; A=0.1; K=1;
beta_1=0.3; beta_2=0.15;
gamma_1=0.02; gamma_2=0.05;
alpha_1=0.02; alpha_2=0.15;
sigma=3;

F = @(x) [
    r*x(1)*(1-x(1)/K)*(x(1)/A-1) - (beta_1*x(2) + beta_2*x(3))*x(1)/(x(1)+x(2)+x(3)) + gamma_1*x(2) + gamma_2*x(3)-D*x(1);
    x(2)*(beta_1*x(1)/(x(1)+x(2)+x(3))- gamma_1 - sigma*beta_2*x(3)/(x(1)+x(2)+x(3)))-D*x(2)-alpha_1*x(2);
    x(3)*(beta_2*x(1)/(x(1)+x(2)+x(3)) - gamma_2 + sigma*beta_2*x(2)/(x(1)+x(2)+x(3)))-D*x(3)-alpha_2*x(3);
];

% initial guess
x0 = [0.8,1,0.2];

options = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');
x_eq = fsolve(F, x0);

S_eq = x_eq(1);
I_eq = x_eq(2);
J_eq = x_eq(3);

disp(['S_eq = ', num2str(S_eq)])
disp(['I_eq = ', num2str(I_eq)])
disp(['J_eq = ', num2str(J_eq)])