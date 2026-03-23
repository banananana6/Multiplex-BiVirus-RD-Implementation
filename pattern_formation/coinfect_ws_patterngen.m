% pattern generation on Watts-Strogatz network

output_folder = ""; % enter output folder here
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

N=1600; p=0.05;

G1 = WattsStrogatz(N, 4, p); % each layer has avg deg of 8
G2 = WattsStrogatz(N, 4, p);
G3 = WattsStrogatz(N, 4, p);
G4 = WattsStrogatz(N, 4, p);

% saving the graph
% ---

graphS = figure('Visible','off', 'Position', [100, 100, 500, 400]);
p1 = plot(G1, 'NodeColor', 'k', ...
    'NodeLabel', {}, ...
    'MarkerSize', 1);
axis off;
savefig(graphS, fullfile(output_folder, 'G_S.fig'));
exportgraphics(graphS, fullfile(output_folder, 'G_S.png'), 'Resolution', 600);
close(graphS);

graphI = figure('Visible','off', 'Position', [100, 100, 500, 400]);
p2 = plot(G2, 'NodeColor', 'k', ...
    'NodeLabel', {}, ...
    'MarkerSize', 1);
axis off;
savefig(graphI, fullfile(output_folder, 'G_I.fig'));
exportgraphics(graphI, fullfile(output_folder, 'G_I.png'), 'Resolution', 600);
close(graphI);

graphJ = figure('Visible','off', 'Position', [100, 100, 500, 400]);
p3 = plot(G3, 'NodeColor', 'k', ...
    'NodeLabel', {}, ...
    'MarkerSize', 1);
axis off;
savefig(graphJ, fullfile(output_folder, 'G_J.fig'));
exportgraphics(graphJ, fullfile(output_folder, 'G_J.png'), 'Resolution', 600);
close(graphJ);

graphC = figure('Visible','off', 'Position', [100, 100, 500, 400]);
p3 = plot(G4, 'NodeColor', 'k', ...
    'NodeLabel', {}, ...
    'MarkerSize', 1);
axis off;
savefig(graphC, fullfile(output_folder, 'G_C.fig'));
exportgraphics(graphC, fullfile(output_folder, 'G_C.png'), 'Resolution', 600);
close(graphC);

%---

A1 = adjacency(G1);
A2 = adjacency(G2);
A3 = adjacency(G3);
A4 = adjacency(G4);

D=0.005; r=0.1; A=0.1; K=1;
beta_1=0.3; beta_2=0.5; beta_10=0.15; beta_02=0.1; beta_12=0.05;
gamma_1=0.02; gamma_2=0.05;
alpha_1=0.08; alpha_2=0.08; alpha_12=0.12;

dS=0.1; dSI=-0.2; dSJ=-0.2; dI=0.1; dJ=0.9; dC=0; dSC=0; dIC=0; dJC=0; dCI=0; dCJ=0;

L1 = A1-diag(sum(A1, 2)); 
L2 = A2-diag(sum(A2, 2)); 
L3 = A3-diag(sum(A3, 2));
L4 = A4-diag(sum(A4, 2));

S_eq = 0.50787
I_eq = 0.29974
J_eq = 0.27653
C_eq = 0.40335

S0 = S_eq+1e-5*randn(N, 1);
I0 = I_eq+1e-5*randn(N, 1);
J0 = J_eq+1e-5*randn(N, 1);
C0 = C_eq+1e-5*randn(N, 1);

y0 = [S0; I0; J0; C0];

num_t = 100;
time_bet = 2;

amplitude_overall = zeros(1, num_t+1);
time_vals = zeros(1, num_t+1);

for i = 1:num_t
    tspan = [0 time_bet];
    
    ode_wrap = @(t, y) coinfect_diff(t, y, L1, L2, L3, L4, r, K, A, D, beta_1, beta_2, beta_10, beta_02, beta_12, gamma_1, gamma_2, alpha_1, alpha_2, alpha_12, dS, dI, dJ, dC, dSI, dSJ, dSC, dIC, dJC, dCI, dCJ, N);
    [t, y] = ode113(ode_wrap, tspan, y0);

    S = y(end, 1:N)';
    I = y(end, N+1:2*N)';
    J = y(end, 2*N+1:3*N)';
    C = y(end, 3*N+1:end)';

    if mod(i, 5) == 0

        cbSize = 14;
        
        figS = figure('Visible','off', 'Position', [100, 100, 550, 400]);
        plot(S, 'b', 'LineWidth', 1.5);
        xlabel('Node Index', 'Interpreter', 'latex', 'FontSize', cbSize);
        ylabel('$S_i$', 'Interpreter', 'latex', 'FontSize', cbSize);
        grid off;
        set(gca, 'FontSize', cbSize, 'TickLabelInterpreter', 'latex');
        exportgraphics(figS, fullfile(output_folder, sprintf('S_t%03d.png', i*time_bet)), 'Resolution', 600);
        close(figS);
        
        figI = figure('Visible','off', 'Position', [100, 100, 550, 400]);
        plot(I, 'r', 'LineWidth', 1.5);
        xlabel('Node Index', 'Interpreter', 'latex', 'FontSize', cbSize);
        ylabel('$I_i$', 'Interpreter', 'latex', 'FontSize', cbSize);
        grid off;
        set(gca, 'FontSize', cbSize, 'TickLabelInterpreter', 'latex');
        exportgraphics(figI, fullfile(output_folder, sprintf('I_t%03d.png', i*time_bet)), 'Resolution', 600);
        close(figI);
        
        figJ = figure('Visible','off', 'Position', [100, 100, 550, 400]);
        plot(J, 'r', 'LineWidth', 1.5);
        xlabel('Node Index', 'Interpreter', 'latex', 'FontSize', cbSize);
        ylabel('$J_i$', 'Interpreter', 'latex', 'FontSize', cbSize);
        grid off;
        set(gca, 'FontSize', cbSize, 'TickLabelInterpreter', 'latex');
        exportgraphics(figJ, fullfile(output_folder, sprintf('J_t%03d.png', i*time_bet)), 'Resolution', 600);
        close(figJ);
        
        figC = figure('Visible','off', 'Position', [100, 100, 550, 400]);
        plot(C, 'b', 'LineWidth', 1.5);
        xlabel('Node Index', 'Interpreter', 'latex', 'FontSize', cbSize);
        ylabel('$C_i$', 'Interpreter', 'latex', 'FontSize', cbSize);
        grid off;
        set(gca, 'FontSize', cbSize, 'TickLabelInterpreter', 'latex');
        exportgraphics(figC, fullfile(output_folder, sprintf('C_t%03d.png', i*time_bet)), 'Resolution', 600);
        close(figC);
    end


    amplitude_overall(i+1) = sqrt(sum((S-S_eq).^2+(I-I_eq).^2+(J-J_eq).^2+(C-C_eq).^2));
    time_vals(i+1) = i*time_bet;

    y0 = [S; I; J; C];
end

close all;
figAmpOverall = figure('Visible', 'on');
clf; hold on;

figAmpOverall = figure;
plot(time_vals, amplitude_overall, '-', 'LineWidth', 1);
xlabel('Time');
ylabel('Amplitude');
grid on;

savefig(figAmpOverall, fullfile(output_folder, 'amp_overall'));

close(figAmpOverall);

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

function h = WattsStrogatz(N, K, beta)
% N is # of nodes, 2*K is each node's deg, beta is prob of rewiring
s = repelem((1:N)', 1, K);
t = s + repmat(1:K, N, 1);
t = mod(t-1, N) + 1;
s = s(:); t = t(:);

for idx=1:length(s)
    if rand<beta
        newTarget = randi(N);
        while newTarget==s(idx) || any(t(s==s(idx))==newTarget)
            newTarget = randi(N);
        end
        t(idx) = newTarget;
    end
end

h = graph(s, t);
end