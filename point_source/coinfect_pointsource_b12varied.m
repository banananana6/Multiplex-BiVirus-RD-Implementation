% point source infections for distinct beta12 values (coinfection)

output_folder = "";
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

rows=40; cols=40; N=rows*cols;

k_1=1; k_2=1; k_3=1; k_4=1;

A1=generate_LA2D(rows, cols, k_1);
A2=generate_LA2D(rows, cols, k_2);
A3=generate_LA2D(rows, cols, k_3);
A4=generate_LA2D(rows, cols, k_4);

D=0.005; r=0.1; A=0.1; K=1;
beta_1=0.3; beta_2=0.4; beta_10=0.2; beta_02=0.3; beta_12=0.05;
gamma_1=0.1; gamma_2=0.05;
alpha_1=0.05; alpha_2=0.15; alpha_12=0.25;

dS=0.5; dSI=0.1; dSJ=0.1; dI=0.4; dJ=0.2; dC=0.1; dSC=0.1; dIC=0; dJC=0; dCI=0; dCJ=0;

L1 = A1-diag(sum(A1, 2)); 
L2 = A2-diag(sum(A2, 2)); 
L3 = A3-diag(sum(A3, 2));
L4 = A4-diag(sum(A4, 2));

all_b12 = [0.01,0.05,0.1,0.14];

disp(all_b12);

colors = lines(length(all_b12));

figI = figure('Visible', 'off'); hold on;
figJ = figure('Visible', 'off'); hold on;
figC = figure('Visible', 'off'); hold on;

num_t=700;
time_bet=1;

for d_idx=1:length(all_b12)
    d_idx
    beta_12=all_b12(d_idx);

    S0=ones(N, 1);
    I0=zeros(N, 1); J0=zeros(N, 1); C0=zeros(N,1);
    I0(1)=0.05; J0(1600)=0.05;

    y0 = [S0;I0;J0;C0];

    nodes_I = zeros(1, num_t+1);
    nodes_J = zeros(1, num_t+1);
    nodes_C = zeros(1, num_t+1);
    time_vals = zeros(1, num_t+1);

    nodes_I(1) = sum(I0 > 0.01)/N;
    nodes_J(1) = sum(J0 > 0.01)/N;
    nodes_C(1) = sum(C0 > 0.005)/N;
    time_vals(1) = 0;

    for i=1:num_t
 
        tspan = [0 time_bet];
        ode_wrap = @(t, y) coinfect_diff(t, y, L1, L2, L3, L4, r, K, A, D, beta_1, beta_2, ...
           beta_10, beta_02, beta_12, gamma_1, gamma_2, ...
           alpha_1, alpha_2, alpha_12, dS, dI, dJ, dC, dSI, dSJ, dSC, dIC, dJC, dCI, dCJ, N);
        [t, y] = ode113(ode_wrap, tspan, y0);

        S = y(end, 1:N)';
        I = y(end, N+1:2*N)';
        J = y(end, 2*N+1:3*N)';
        C = y(end, 3*N+1:end)';

        nodes_I(i+1)=sum(I > 0.01)/N;
        nodes_J(i+1)=sum(J > 0.01)/N;
        nodes_C(i+1)=sum(C > 0.005)/N;

        time_vals(i+1)=i*time_bet;

        y0 = [S;I;J;C];
    end

    figure(figI); 
    plot(time_vals, nodes_I, '-', 'Color', colors(d_idx,:), 'LineWidth', 1);

    figure(figJ); 
    plot(time_vals, nodes_J, '-', 'Color', colors(d_idx,:), 'LineWidth', 1);

    figure(figC); 
    plot(time_vals, nodes_C, '-', 'Color', colors(d_idx,:), 'LineWidth', 1);
end

figure(figI); xlabel('Time'); ylabel('I-Spread Index');
legend(compose('b_{12} = %.3f', all_b12), 'Location', 'best'); grid on;
saveas(figI, fullfile(output_folder, 'b12_I'));
close(figI);

figure(figJ); xlabel('Time'); ylabel('J-Spread Index');
legend(compose('b_{12} = %.3f', all_b12), 'Location', 'best'); grid on;
saveas(figJ, fullfile(output_folder, 'b12_J'));
close(figJ);

figure(figC); xlabel('Time'); ylabel('Co-Spread Index');
legend(compose('b_{12} = %.3f', all_b12), 'Location', 'best'); grid on;
saveas(figC, fullfile(output_folder, 'b12_C'));
close(figC);

function dydt = coinfect_diff(t, y, L1, L2, L3, L4, r, K, A, D, beta_1, beta_2, beta_10, beta_02, beta_12, gamma_1, gamma_2, alpha_1, alpha_2, alpha_12, dS, dI, dJ, dC, dSI, dSJ, dSC, dIC, dJC, dCI, dCJ, N)
    S = y(1:N);
    I = y(N+1:2*N);
    J = y(2*N+1:3*N);
    C = y(3*N+1:end);

    f = r*S.*(1-S/K).*(S/A-1)-(beta_1*I + beta_2*J+(beta_10+beta_02+beta_12)*C).*S./(S+I+J+C) ...  
    +gamma_1*I+gamma_2*J-D*S;
    g = (beta_1*I+beta_10*C).*S./(S+I+J+C)-(beta_2*J+beta_02*C+beta_12*C).*I./(S+I+J+C)+(gamma_2*C-gamma_1*I)-D*I-alpha_1*I;
    h=(beta_2*J+beta_02*C).*S./(S+I+J+C)-(beta_1*I+beta_10*C+beta_12*C).*J./(S+I+J+C)+(gamma_1*C-gamma_2*J)-D*J-alpha_2*J;
    l=(beta_12*C).*S./(S+I+J+C)+(beta_2*J+beta_02*C+beta_12*C).*I./(S+I+J+C)+(beta_1*I+beta_10*C+beta_12*C).*J./(S+I+J+C)...
        -(gamma_1+gamma_2)*C-D*C-alpha_12*C;

    dSdt = f+dS*(L1*S)+dSI*(L2*I)+dSJ*(L3*J)+dSC*(L4*C);
    dIdt = g+dI*(L2*I)+dIC*(L4*C);
    dJdt = h+dJ*(L3*J)+dJC*(L4*C);
    dCdt = l+dC*(L4*C)+dCI*(L2*I)+dCJ*(L3*J);

    dydt = [dSdt; dIdt; dJdt; dCdt];
end

function A = generate_LA2D(rows, cols, range)
    N = rows*cols;
    A = zeros(N);

    for idx=1:N
        yi = ceil(idx/cols);
        xi = idx-(yi-1)*cols;

        neighbors = [];
        for dx = -range:range
            for dy = -range:range
                if dx==0 && dy==0
                    continue;
                end
                d2 = dx^2+dy^2;
                if d2<=range^2 && d2~=8
                    xj = xi+dx;
                    yj = yi+dy;
                    if xj>=1 && xj<=cols && yj>=1 && yj<=rows
                        jdx = (yj-1)*cols+xj;
                        neighbors(end+1) = jdx;
                    end
                end
            end
        end

        neighbors = sort(neighbors);
        A(idx, neighbors) = 1;
    end
end