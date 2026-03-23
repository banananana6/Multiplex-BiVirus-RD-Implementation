% point source spread for distinct outbreak timing (superinfection)

output_folder = "";
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

rows=40; cols=40; N=rows*cols;
k_1=1; k_2=1; k_3=1;

A1 = generate_LA2D(rows, cols, k_1);
A2 = generate_LA2D(rows, cols, k_2);
A3 = generate_LA2D(rows, cols, k_3);

D=0.005; r=0.1; A=0.1; K=1;
beta_1=0.5; beta_2=0.4;
gamma_1=0.2; gamma_2=0.1;
alpha_1=0.01; alpha_2=0.05; sigma=0.9;

dS=0.3; dSI=0.1; dSJ=0.1; dI=0.3; dJ=0.1;

L1 = A1-diag(sum(A1, 2)); 
L2 = A2-diag(sum(A2, 2)); 
L3 = A3-diag(sum(A3, 2));

all_times=linspace(5, 80, 16);

disp(all_times);

colors=lines(length(all_times));

num_t=500;
time_bet=1;

I_peaks = zeros(1, length(all_times));
I_peak_times = zeros(1, length(all_times));
J_peaks = zeros(1, length(all_times));
J_peak_times = zeros(1, length(all_times));

for d_idx=1:length(all_times)
    d_idx
    d=1600;

    intro_time = all_times(d_idx);

    S0=ones(N, 1);
    I0=zeros(N, 1); J0=zeros(N, 1);
    I0(1)=0.05;

    y0 = [S0;I0;J0];

    nodes_I = zeros(1, num_t+1);
    nodes_J = zeros(1, num_t+1);
    time_vals = zeros(1, num_t+1);

    nodes_I(1) = sum(I0>0.01)/N;
    nodes_J(1) = sum(J0>0.01)/N;
    time_vals(1) = 0;

    sat_time = NaN;

    for i = 1:num_t
 
        tspan = [0 time_bet];
        ode_wrap = @(t, y) superinfect_diff(t, y, L1, L2, L3, r, K, A, D, beta_1, beta_2, ...
            gamma_1, gamma_2, alpha_1, alpha_2, sigma, dS, dSI, dSJ, dI, dJ, N);
        [t, y] = ode113(ode_wrap, tspan, y0);

        S=y(end, 1:N)'; I=y(end, N+1:2*N)'; J=y(end, 2*N+1:end)';

        if i==intro_time
             J(1600)=0.05;
        end

        nodes_I(i+1) = sum(I>0.01)/N;
        nodes_J(i+1) = sum(J>0.01)/N;
        time_vals(i+1) = i*time_bet;

        y0 = [S; I; J];
    end

    [I_peaks(d_idx), idx_peak]=max(nodes_I);
    I_peak_times(d_idx)=time_vals(idx_peak);

    [J_peaks(d_idx), idx_peak]=max(nodes_J);
    J_peak_times(d_idx)=time_vals(idx_peak);
end

figure;
scatter(all_times, I_peaks, 40, 'filled');
xlabel('t_d');
ylabel('I_1 Saturation Time');
filename = sprintf('I1_peak_vs_td');
savefig(fullfile(output_folder, filename));

figure;
scatter(all_times, I_peak_times, 40, 'filled');
xlabel('t_d');
ylabel('I_1 Saturation Time');
filename = sprintf('I1_time_vs_td');
savefig(fullfile(output_folder, filename));

figure;
scatter(all_times, J_peaks, 40, 'filled');
xlabel('t_d');
ylabel('I_2 Saturation Time');
filename = sprintf('I2_peak_vs_td');
savefig(fullfile(output_folder, filename));

figure;
scatter(all_times, J_peak_times, 40, 'filled');
xlabel('t_d');
ylabel('I_2 Saturation Time');
filename = sprintf('I2_time_vs_td');
savefig(fullfile(output_folder, filename));

function dydt = superinfect_diff(t,y,L1,L2,L3,r,K,A,D,beta_1,beta_2,gamma_1,gamma_2,alpha_1,alpha_2,sigma,dS,dSI,dSJ,dI,dJ,N)
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