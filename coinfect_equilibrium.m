D=0.005; r=0.1; A=0.1; K=1;
beta_1=0.3; beta_2=0.5; beta_10=0.15; beta_02=0.1; beta_12=0.05;
gamma_1=0.02; gamma_2=0.05;
alpha_1=0.08; alpha_2=0.08; alpha_12=0.12;

F = @(x) [
    r*x(1)*(1-x(1)/K)*(x(1)/A-1) - (beta_1*x(2) + beta_2*x(3) + (beta_10+beta_02+beta_12)*x(4))*x(1)/(x(1)+x(2)+x(3)+x(4)) ...
    + gamma_1*x(2) + gamma_2*x(3) - D*x(1);
    
    (beta_1*x(2) + beta_10*x(4))*x(1)/(x(1)+x(2)+x(3)+x(4)) - (beta_2*x(3) + beta_02*x(4) + beta_12*x(4))*x(2)/(x(1)+x(2)+x(3)+x(4)) ...
    + (gamma_2*x(4) - gamma_1*x(2)) - D*x(2) - alpha_1*x(2);

    (beta_2*x(3) + beta_02*x(4))*x(1)/(x(1)+x(2)+x(3)+x(4)) - (beta_1*x(2) + beta_10*x(4) + beta_12*x(4))*x(3)/(x(1)+x(2)+x(3)+x(4)) ...
    + (gamma_1*x(4) - gamma_2*x(3)) - D*x(3) - alpha_2*x(3);
    
    (beta_12*x(4))*x(1)/(x(1)+x(2)+x(3)+x(4)) + (beta_2*x(3) + beta_02*x(4) + beta_12*x(4))*x(2)/(x(1)+x(2)+x(3)+x(4)) ...
    + (beta_1*x(2) + beta_10*x(4) + beta_12*x(4))*x(3)/(x(1)+x(2)+x(3)+x(4)) ...
    - (gamma_1 + gamma_2)*x(4) - D*x(4) - alpha_12*x(4);
];

x0 = [0.5, 0.2, 0.2, 0.2];

options = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');
x_eq = fsolve(F, x0, options);

S_eq = x_eq(1);
I_eq = x_eq(2);
J_eq = x_eq(3);
C_eq = x_eq(4);

disp(['S_eq = ', num2str(S_eq)])
disp(['I_eq = ', num2str(I_eq)])
disp(['J_eq = ', num2str(J_eq)])
disp(['C_eq = ', num2str(C_eq)])