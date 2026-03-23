% superinfection parameter sensitivity when sigma=2

output_folder = "";
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

rows=40; cols=40; N=rows*cols;
k_1=2; k_2=2; k_3=2;

A1 = generate_LA2D(rows, cols, k_1);
A2 = generate_LA2D(rows, cols, k_2);
A3 = generate_LA2D(rows, cols, k_3);

D=0.005; r=0.1; A=0.1; K=1;
beta_1=0.3; beta_2=0.15;
gamma_1=0.02; gamma_2=0.05;
alpha_1=0.02; alpha_2=0.15;
sigma=2;

S_eq = 0.27646
I_eq = 0.469
J_eq = 0.14317

L1 = A1-diag(sum(A1, 2)); 
L2 = A2-diag(sum(A2, 2)); 
L3 = A3-diag(sum(A3, 2));

dS=0.1; dI=0.01; dJ=4.8;

dSI_vals = -0.3:0.05:0.3;
dSJ_vals = -0.3:0.05:0.3;

tspan = [0 50];
for dSI=dSI_vals
    for dSJ=dSJ_vals
        S0 = S_eq+1e-3*randn(N, 1);
        I0 = I_eq+1e-3*randn(N, 1);
        J0 = J_eq+1e-3*randn(N, 1);
        y0 = [S0; I0; J0];

        ode_wrap = @(t, y) superinfect_diff(t, y, L1, L2, L3, ...
            r, K, A, D, beta_1, beta_2, gamma_1, gamma_2, ...
            alpha_1, alpha_2, sigma, ...
            dS, dSI, dSJ, dI, dJ, N);
        [t, y] = ode113(ode_wrap, tspan, y0);

        S_final = y(end, 1:N)';
        I_final = y(end, N+1:2*N)';
        J_final = y(end, 2*N+1:end)';

        S_plot = reshape(S_final, [cols, rows])';
        I_plot = reshape(I_final, [cols, rows])';
        J_plot = reshape(J_final, [cols, rows])';

        fig = figure('Visible', 'off');
        imagesc(S_plot); colormap('parula'); colorbar;
        axis equal tight; title('S');
        set(gca, 'FontSize', 12);
        saveas(fig, fullfile(output_folder, sprintf('S_dSI_%.2f_dSJ_%.2f.png', dSI, dSJ)));
        close(fig);

        fig = figure('Visible', 'off');
        imagesc(I_plot); colormap('hot'); colorbar;
        axis equal tight; title('I');
        set(gca, 'FontSize', 12);
        saveas(fig, fullfile(output_folder, sprintf('I_dSI_%.2f_dSJ_%.2f.png', dSI, dSJ)));
        close(fig);

        fig = figure('Visible', 'off');
        imagesc(J_plot); colormap('hot'); colorbar;
        axis equal tight; title('J');
        set(gca, 'FontSize', 12);
        saveas(fig, fullfile(output_folder, sprintf('J_dSI_%.2f_dSJ_%.2f.png', dSI, dSJ)));
        close(fig);
    end
end


function dydt = superinfect_diff(t, y, L1, L2, L3, r, K, A, D, beta_1, beta_2, gamma_1, gamma_2, alpha_1, alpha_2, sigma, dS, dSI, dSJ, dI, dJ, N)
    S = y(1:N);
    I = y(N+1:2*N);
    J = y(2*N+1:end);

    f = r*S.*(1-S/K).*(S/A-1)-(beta_1*I+beta_2*J).*S./(S+I+J)+gamma_1*I+gamma_2*J-D*S;
    g = I.*(beta_1*S./(S+I+J)-alpha_1-gamma_1-sigma*beta_2*J./(S+I+J))-D*I;
    h = J.*(beta_2*S./(S+I+J)-alpha_2-gamma_2+sigma*beta_2*I./(S+I+J))-D*J;

    dSdt = f+dS*(L1*S)+dSI*(L2*I)+dSJ*(L3*J);
    dIdt = g+dI*(L2*I);
    dJdt = h+dJ*(L3*J);

    dydt = [dSdt; dIdt; dJdt];
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