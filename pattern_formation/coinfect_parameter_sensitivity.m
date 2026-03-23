output_folder = ""; % output folder here
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

rows=40; cols=40; N=rows*cols;

k_1=1; k_2=1; k_3=1; k_4=1; % LA4-LA4-LA4-LA4

A1 = generate_LA2D(rows, cols, k_1);
A2 = generate_LA2D(rows, cols, k_2);
A3 = generate_LA2D(rows, cols, k_3);
A4 = generate_LA2D(rows, cols, k_4);

D=0.005; r=0.1; A=0.1; K=1;
beta_1=0.3; beta_2=0.5; beta_10=0.15; beta_02=0.1; beta_12=0.05;
gamma_1=0.02; gamma_2=0.05;
alpha_1=0.08; alpha_2=0.08; alpha_12=0.12;

dS=0.1; dSI=-0.2; dSJ=-0.2; dSC=-0.2;
dI=0.1; dJ=0.08; 
dCI=0; dCJ=0; dC=0.5;

dIC_vals=-0.3:0.05:0.3;
dJC_vals=-0.3:0.05:0.3;

S_eq = 0.50787
I_eq = 0.29974
J_eq = 0.27653
C_eq = 0.40335

L1 = A1-diag(sum(A1, 2)); 
L2 = A2-diag(sum(A2, 2)); 
L3 = A3-diag(sum(A3, 2));
L4 = A4-diag(sum(A4, 2));

tspan = [0 110];
for dIC=dIC_vals
    for dJC=dJC_vals
        S0 = S_eq+1e-5*randn(N, 1);
        I0 = I_eq+1e-5*randn(N, 1);
        J0 = J_eq+1e-5*randn(N, 1);
        C0 = C_eq+1e-5*randn(N, 1);
        y0 = [S0; I0; J0; C0];

        ode_wrap = @(t, y) coinfect_diff(t, y, L1, L2, L3, L4, r, K, A, D, beta_1, beta_2, beta_10, beta_02, beta_12, gamma_1, gamma_2, alpha_1, alpha_2, alpha_12, dS, dI, dJ, dC, dSI, dSJ, dSC, dIC, dJC, dCI, dCJ, N);
        [t, y] = ode113(ode_wrap, tspan, y0);

        S_final = y(end, 1:N)';
        I_final = y(end, N+1:2*N)';
        J_final = y(end, 2*N+1:3*N)';
        C_final = y(end, 3*N+1:end);

        S_plot = reshape(S_final, [cols, rows])';
        I_plot = reshape(I_final, [cols, rows])';
        J_plot = reshape(J_final, [cols, rows])';
        C_plot = reshape(C_final,[cols,rows])';

        figS = figure('Visible','off');
        imagesc(S_plot);
        colormap(figS,'turbo');
        c = colorbar; 
        c.TickLabelInterpreter = 'latex';
        c.FontSize = cbSize;
        axis equal tight;
        axis off;
        exportgraphics(figS, fullfile(output_folder, sprintf('S_dIC_%.2f_dJC_%.2f_d44_%.2f.png', dIC, dJC, dC)), 'Resolution', 600);
        close(figS);

        figI = figure('Visible','off');
        imagesc(I_plot);
        colormap(figI,'turbo');
        c = colorbar; 
        c.TickLabelInterpreter = 'latex';
        c.FontSize = cbSize;
        axis equal tight;
        axis off;
        exportgraphics(figI, fullfile(output_folder, sprintf('I_dIC_%.2f_dJC_%.2f_dC_%.2f.png', dIC, dJC, dC)), 'Resolution', 600);
        close(figI);
            
        figJ = figure('Visible','off');
        imagesc(J_plot);
        colormap(figJ,'turbo');
        c = colorbar; 
        c.TickLabelInterpreter = 'latex';
        c.FontSize = cbSize;
        axis equal tight;
        axis off;
        exportgraphics(figJ, fullfile(output_folder, sprintf('J_dIC_%.2f_dJC_%.2f_dC_%.2f.png', dIC, dJC, dC)), 'Resolution', 600);
        close(figJ);
    
        figC = figure('Visible','off');
        imagesc(C_plot);
        colormap(figC,'turbo');
        c = colorbar; 
        c.TickLabelInterpreter = 'latex';
        c.FontSize = cbSize;
        axis equal tight;
        axis off;
        exportgraphics(figC, fullfile(output_folder, sprintf('C_dIC_%.2f_dJC_%.2f_dC_%.2f.png', dIC, dJC,d C)), 'Resolution', 600);
        close(figC);
    end
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