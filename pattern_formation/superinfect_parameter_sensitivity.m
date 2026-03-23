% superinfect parameter sensitivity

rows=40; cols=40; N=rows*cols;
k_1=2; k_2=2; k_3=2;

A1 = generate_LA2D(rows, cols, k_1);
A2 = generate_LA2D(rows, cols, k_2);
A3 = generate_LA2D(rows, cols, k_3);

D=0.005; r=0.1; A=0.1; K=1;
beta_1=0.3; beta_2=0.15;
gamma_1=0.02; gamma_2=0.05;
alpha_1=0.02; alpha_2=0.15;
sigma=3;

S_eq = 0.88318
I_eq = 0.53803
J_eq = 0.40606

L1 = A1-diag(sum(A1, 2)); 
L2 = A2-diag(sum(A2, 2)); 
L3 = A3-diag(sum(A3, 2));

S0 = S_eq + 1e-3*randn(N, 1);
I0 = I_eq + 1e-3*randn(N, 1);
J0 = J_eq + 1e-3*randn(N, 1);
y0 = [S0; I0; J0];

d11_vals = 0.1:1:0.1;
d12_vals = -0.3:0.05:0.3;
d13_vals = -0.3:0.05:0.3;
d22_vals = 0.01:1:0.01;
d33_vals = 4.8:1:4.8;

tspan = [0 200];
for d11 = d11_vals
    for d12 = d12_vals
        for d13 = d13_vals
            for d22 = d22_vals
                for d33 = d33_vals
                    S0 = S_eq + 1e-3 * randn(N, 1);
                    I0 = I_eq + 1e-3 * randn(N, 1);
                    J0 = J_eq + 1e-3 * randn(N, 1);
                    y0 = [S0; I0; J0];

                    ode_wrap = @(t, y) superinfect_diff(t, y, L1, L2, L3, ...
                        r, K, A, D, beta_1, beta_2, gamma_1, gamma_2, ...
                        alpha_1, alpha_2, sigma, ...
                        d11, d12, d13, d22, d33, N);
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
                    saveas(fig, fullfile(output_folder, sprintf('S_d12_%.2f_d13_%.2f.png',d12, d13)));
                    close(fig);

                    fig = figure('Visible', 'off');
                    imagesc(I_plot); colormap('hot'); colorbar;
                    axis equal tight; title('I');
                    set(gca, 'FontSize', 12);
                    saveas(fig, fullfile(output_folder, sprintf('I_d12_%.2f_d13_%.2f.png',d12, d13)));
                    close(fig);

                    fig = figure('Visible', 'off');
                    imagesc(J_plot); colormap('hot'); colorbar;
                    axis equal tight; title('J');
                    set(gca, 'FontSize', 12);
                    saveas(fig, fullfile(output_folder, sprintf('J_d12_%.2f_d13_%.2f.png',d12, d13)));
                    close(fig);
                end
            end
        end
    end
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