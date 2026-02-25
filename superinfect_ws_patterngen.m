N=1600;
p = 0.05;

output_folder = ""; % output folder here
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

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

D=0.005; r=0.1; A=0.1; K=1;
beta_1=0.3; beta_2=0.15;
gamma_1=0.02; gamma_2=0.05;
alpha_1=0.02; alpha_2=0.15;
sigma=3;

d11=0.1; d12=-0.2; d13=-0.2; d22=0.01; d33=4;

S_eq = 0.88318
I_eq = 0.53803
J_eq = 0.40606

G1 = WattsStrogatz(N, 6, p);
G2 = G1;
G3 = G1;

A1 = adjacency(G1);
A2 = adjacency(G2);
A3 = adjacency(G3);

L1 = A1-diag(sum(A1, 2)); 
L2 = A2-diag(sum(A2, 2)); 
L3 = A3-diag(sum(A3, 2));

S0 = S_eq + 1e-3*randn(N, 1);
I0 = I_eq + 1e-3*randn(N, 1);
J0 = J_eq + 1e-3*randn(N, 1);
y0 = [S0; I0; J0];

num_t=1000;
time_bet=2;

amplitude_overall = zeros(1, num_t+1);
time_vals = zeros(1, num_t+1);

for i = 1:num_t
    tspan = [0 time_bet];
    
    ode_wrap = @(t, y) superinfect_diff(t, y, L1, L2, L3, r,K,A,D,beta_1, beta_2, gamma_1, gamma_2, alpha_1, alpha_2, sigma, d11, d12, d13, d22, d33, N);
    [t, y] = ode113(ode_wrap, tspan, y0);

    S = y(end, 1:N)';
    I = y(end, N+1:2*N)';
    J = y(end, 2*N+1:end)';

    if mod(i, 10) == 0
        figS = figure('Visible', 'off');
        set(figS, 'Units', 'inches', 'Position', [1 1 5 4]); % width=5, height=4
        scatter(1:numel(S), S, 10, 'b', 'filled');
        xlabel('Node Index','Interpreter','latex');
        ylabel('$S_i$','Interpreter','latex');
        set(gca, 'FontSize', 12, 'TickLabelInterpreter','latex', ...
                 'LineWidth', 1, 'Box','on');
        exportgraphics(figS, fullfile(output_folder, sprintf('S_t%03d.png', i*time_bet)), ...
                       'Resolution', 600);
        close(figS);
        
        figI = figure('Visible', 'off');
        set(figI, 'Units', 'inches', 'Position', [1 1 5 4]); % width=5, height=4
        scatter(1:numel(I), I, 10, 'r', 'filled');
        xlabel('Node Index','Interpreter','latex');
        ylabel('$I_i$','Interpreter','latex');
        set(gca, 'FontSize', 12, 'TickLabelInterpreter','latex', ...
                 'LineWidth', 1, 'Box','on');
        exportgraphics(figI, fullfile(output_folder, sprintf('I_t%03d.png', i*time_bet)), ...
                       'Resolution', 600);
        close(figI);
        
        figJ = figure('Visible', 'off');
        set(figJ, 'Units', 'inches', 'Position', [1 1 5 4]); % width=5, height=4
        scatter(1:numel(J), J, 25, 'g', 'filled');
        xlabel('Node Index','Interpreter','latex');
        ylabel('$J_i$','Interpreter','latex');
        set(gca, 'FontSize', 12, 'TickLabelInterpreter','latex', ...
                 'LineWidth', 1, 'Box','on');
        exportgraphics(figJ, fullfile(output_folder, sprintf('J_t%03d.png', i*time_bet)), ...
                       'Resolution', 600);
        close(figJ);
    end

    amplitude_overall(i+1) = sqrt(sum((S - S_eq).^2 + (I - I_eq).^2 + (J - J_eq).^2));
    time_vals(i+1) = i*time_bet;

    y0 = [S; I; J];
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

function dydt = superinfect_diff(t, y, L1, L2, L3, r, K, A, D, beta_1, beta_2, gamma_1, gamma_2, alpha_1, alpha_2, sigma, d11, d12, d13, d22, d33, N)
    S = y(1:N);
    I = y(N+1:2*N);
    J = y(2*N+1:end);

    f = r*S.*(1-S/K).*(S/A-1)-(beta_1*I+beta_2*J).*S./(S+I+J)+gamma_1*I+gamma_2*J-D*S;
    g = I.*(beta_1*S./(S+I+J)-gamma_1-sigma*beta_2*J./(S+I+J))-D*I-alpha_1*I;
    h = J.*(beta_2*S./(S+I+J)-gamma_2+sigma*beta_2*I./(S+I+J))-D*J-alpha_2*J;

    dSdt = f+d11*(L1*S)+d12*(L2*I)+d13*(L3*J);
    dIdt = g+d22*(L2*I);
    dJdt = h+d33*(L3*J);

    dydt = [dSdt; dIdt; dJdt];
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