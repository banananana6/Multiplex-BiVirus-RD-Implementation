% finding beta_12 threshold for each alpha_12 value in a range
% (coinfection)

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

params.L1=L1; params.L2=L2; params.L3=L3; params.L4=L4;
params.r=r; params.K=K; params.A=A; params.D=D;
params.beta_1=beta_1; params.beta_2=beta_2;
params.beta_10=beta_10; params.beta_02=beta_02;
params.gamma_1=gamma_1; params.gamma_2=gamma_2;
params.alpha_1=alpha_1; params.alpha_2=alpha_2;
params.dS=dS; params.dSI=dSI; params.dSJ=dSJ; params.dSC=dSC;
params.dI=dI; params.dJ=dJ; params.dC=dC; params.dIC=dIC;
params.dJC=dJC; params.dCI=dCI; params.dCJ=dCJ;
params.N=N;
params.num_t=1000;
params.time_bet=1;

S0 = ones(N, 1); I0 = zeros(N, 1); J0 = zeros(N, 1); C0 = zeros(N, 1);
I0(1)=0.05; J0(1600)=0.05;
params.y0 = [S0; I0; J0; C0];

alpha12_vals = 0.36:0.02:0.66;  
b12_thresholds = zeros(size(alpha12_vals));

for idx = 1:length(alpha12_vals)
    fprintf('Computing for alpha_12 = %.3f...\n', alpha12_vals(idx));
    params.alpha_12=alpha12_vals(idx);   
    b12_thresholds(idx)=find_b12_threshold([0, 0.5], 1e-3, 20, params);
end

figure;
scatter(alpha12_vals, b12_thresholds, 40, 'filled');
xlabel('\alpha_{12}');
ylabel('b_{12} Threshold');
filename = sprintf('b12thres_vs_alpha_new');
savefig(fullfile(output_folder, filename));

function b12_thresh = find_b12_threshold(b12_range, tol, max_iter, params)
    low=b12_range(1);
    high=b12_range(2);
    b12_thresh=low;

    for iter=1:max_iter
        mid=(low+high)/2;
        [C_peak, C_final] = simulate_coinfect_C(mid, params);
        fprintf('Iter %d: b12=%.4f, C_peak=%.3f, C_final=%.3f\n', iter, mid, C_peak, C_final);

        if C_final==0
            b12_thresh=mid;
            low=mid;
        else
            high=mid;
        end

        if abs(high-low)<tol
            break;
        end
    end
end

function [C_peak, C_final] = simulate_coinfect_C(beta_12, params)
    % Run the coinfection model and return C peak and final values

    y0=params.y0;
    num_t=params.num_t;
    time_bet=params.time_bet;

    nodes_C=zeros(1, num_t+1);
    nodes_C(1)=sum(y0(3*params.N+1:end)>0.005) / params.N;

    for i=1:num_t
        tspan = [0 time_bet];
        ode_wrap = @(t, y) coinfect_diff(t, y, params.L1, params.L2, params.L3, ...
            params.L4, params.r, params.K, params.A, params.D, params.beta_1, params.beta_2, ...
            params.beta_10, params.beta_02, beta_12, params.gamma_1, params.gamma_2, ...
            params.alpha_1, params.alpha_2, params.alpha_12, params.dS, params.dI, params.dJ, ...
            params.dC, params.dSI, params.dSJ, params.dSC, params.dIC, params.dJC, params.dCI, params.dCJ, params.N);
        [~, y] = ode113(ode_wrap, tspan, y0);

        C = y(end, 3*params.N+1:end)';
        C = max(C, 0);
        nodes_C(i+1) = sum(C>0.005)/params.N;

        y0 = y(end,:)';
    end

    C_peak=max(nodes_C);
    C_final=nodes_C(end);
end

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