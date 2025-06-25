clc;
clear;

%% Parameters
N = 1024;             % Code length
n = log2(N);          
epsilon = 0.5;        % BEC erasure probability
EbN0_dB = 2;          % AWGN Eb/N0 in dB
R = 0.5;              % Code rate
EbN0 = 10^(EbN0_dB/10);
sigma = sqrt(1 / (2 * R * EbN0));

%% Bhattacharyya Initialization
Z_bec = zeros(1, N);

% Base initialization
Z_bec(1) = epsilon;

%% Recursive Bhattacharyya Update
for lev = 1:n
    len = 2^lev;
    for i = len/2+1:len
        Z_bec(i)  = 2 * Z_bec(i - len/2) - Z_bec(i - len/2)^2;
    end
    for i = len/2+1:len
        Z_bec(i - len/2)  = Z_bec(i - len/2)^2;
    end
end

%% Final Values
Z_bec_final = Z_bec;

% Sorting
Z_bec_sorted = sort(Z_bec_final);

%% Plotting
figure('Name','Bhattacharyya Parameter Plots (N = 1024)','NumberTitle','off');

subplot(1,2,1);
plot(1:N, Z_bec_final, '.', 'Color', [0 0.45 0.74]);
title('BEC Channel: Unsorted Bhattacharyya');
xlabel('Channel Index');
ylabel('Z(W)');
grid on;

subplot(1,2,2);
plot(1:N, Z_bec_sorted, '.', 'Color', [0.85 0.33 0.1]);
title('BEC Channel: Sorted Bhattacharyya');
xlabel('Sorted Index');
ylabel('Z(W)');
grid on;


