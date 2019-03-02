% Plot the time signatures of each run
clear;clc;close all;

%% First experiment

data_1 = [
    97,49,2883;
    618,314,6684;
    2107,1037,9340;
    5506,2508,11217;
    14464,5153,15612;
    28768,9374,18503;
    51514,15943,19796;
    83141,27007,23126
    ];
N_1 = [101,201,301,401,501,601,701,801];

figure;
for i=1:3
    plot(N_1,data_1(:,i),'-^');
    hold on;
end
grid on;
grid minor;
title('Runtime (ms) against discretisations (Nx,Ny,Nt)');
xlabel('Number of discretisations along Nx, Ny, Nt');
ylabel('Time (ms)');
legend('Serial', 'Parallel (2x1)', 'Parallel (10x10)');
hold off;

%% Second experiment

num_proc = 1:1:10;
time = [3681, 761, 364, 202, 128, 102, 93.6, 85.2, 101.5, 135.05];

figure;
plot(num_proc, time, '-o');
title('Runtime (s) against number of processors along one direction (Px, Py)');
xlabel('Number of processors spanning each direction (Px, Py)');
ylabel('Time (s)');
grid on;
grid minor;