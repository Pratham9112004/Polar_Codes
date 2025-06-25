% Polar Code Encoding and Successive Cancellation Decoding

clearvars;

load('Reliability_Sequence.mat');
blockLength = 512;
n = log2(blockLength);
msgLength = 260;

reliableSet = Reliability_Sequence(Reliability_Sequence <= blockLength);
frozenSet = reliableSet(1:blockLength - msgLength);

Nsim = 250;
EbNo_dB_List = 0:0.5:10;

successRate = zeros(1, length(EbNo_dB_List));
bitErrorRate = zeros(1, length(EbNo_dB_List));
blockErrorRate = zeros(1, length(EbNo_dB_List));
idx = 1;

for Eb_dB = EbNo_dB_List
    successfulBlocks = 0;
    bitErrorsTotal = 0;
    blockErrorsTotal = 0;

    for t = 1:Nsim
        inputBits = randi([0 1], 1, msgLength);
        fullMessage = zeros(1, blockLength);
        fullMessage(reliableSet(blockLength - msgLength + 1:end)) = inputBits;

        encodedBits = polarEncode(fullMessage, n);
        noisyChannelOutput = addBPSKandAWGN(encodedBits, Eb_dB, msgLength, blockLength);

        global L;
        global u_hat;
        global node_counter;

        L = zeros(n + 1, blockLength);
        L(1, :) = noisyChannelOutput;
        u_hat = zeros(n + 1, blockLength);
        node_counter = zeros(n + 1, 1);

        polarDecoder(noisyChannelOutput, 1, n, frozenSet);
        estimatedBits = u_hat(n + 1, reliableSet(blockLength - msgLength + 1:end));

        errors = sum(inputBits ~= estimatedBits);
        if errors == 0
            successfulBlocks = successfulBlocks + 1;
        else
            bitErrorsTotal = bitErrorsTotal + errors;
            blockErrorsTotal = blockErrorsTotal + 1;
        end
    end

    successRate(idx) = successfulBlocks / Nsim;
    bitErrorRate(idx) = bitErrorsTotal / (msgLength * Nsim);
    blockErrorRate(idx) = blockErrorsTotal / Nsim;
    idx = idx + 1;
end

function result = polarEncode(input, stage)
    baseMatrix = [1 0; 1 1];
    G = baseMatrix;
    for i = 1:stage-1
        G = kron(G, baseMatrix);
    end
    result = mod(input * G, 2);
end

function y = addBPSKandAWGN(x, Eb_dB, K, N)
    bpsk = 1 - 2 * x;
    Eb = 10^(Eb_dB / 10);
    Es = (K / N) * Eb;
    noiseStd = sqrt(1 / Es);
    noise = noiseStd * randn(1, N);
    y = bpsk + noise;
end

function polarDecoder(rx, level, depth, frozen)
    global L;
    global u_hat;
    global node_counter;

    f_node = @(a,b) sign(a).*sign(b).*min(abs(a), abs(b));
    g_node = @(a,b,c) b + (1 - 2 * c) .* a;

    idx = node_counter(level);

    if level == depth + 1
        node_counter(level) = node_counter(level) + 1;
        if ismember(idx + 1, frozen)
            u_hat(depth + 1, idx + 1) = 0;
        else
            u_hat(depth + 1, idx + 1) = double(L(depth + 1, idx + 1) < 0);
        end
        return;
    end

    segmentLength = 2^(depth - level + 1);
    a = rx(1:segmentLength/2);
    b = rx(segmentLength/2 + 1:end);
    L(level + 1, segmentLength/2 * idx + 1 : segmentLength/2 * (idx + 1)) = f_node(a, b);
    node_counter(level) = node_counter(level) + 1;

    polarDecoder(f_node(a, b), level + 1, depth, frozen);

    left_bits = u_hat(level + 1, segmentLength/2 * idx + 1 : segmentLength/2 * (idx + 1));

    currentIdx = node_counter(level);
    L(level + 1, segmentLength/2 * currentIdx + 1 : segmentLength/2 * (currentIdx + 1)) = g_node(a, b, left_bits);
    node_counter(level) = node_counter(level) + 1;

    polarDecoder(g_node(a, b, left_bits), level + 1, depth, frozen);
    right_bits = u_hat(level + 1, segmentLength/2 * currentIdx + 1 : segmentLength/2 * (currentIdx + 1));
    left_bits_prev = u_hat(level + 1, segmentLength/2 * idx + 1 : segmentLength/2 * (idx + 1));

    u_hat(level, segmentLength/2 * idx + 1 : segmentLength/2 * (currentIdx + 1)) = ...
        [mod(left_bits_prev + right_bits, 2), right_bits];
end

figure;
plot(EbNo_dB_List, successRate, 'LineWidth', 2);
title("Eb/N0 vs Decoding Success (N,K) = (" + num2str(blockLength) + "," + num2str(msgLength) + ")");
xlabel('Eb/N0 (dB)');
ylabel('Decoding Success Rate');
grid on;

figure;
plot(EbNo_dB_List, bitErrorRate, 'LineWidth', 2);
title("Eb/N0 vs Bit Error Rate (N,K) = (" + num2str(blockLength) + "," + num2str(msgLength) + ")");
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Rate');
grid on;

figure;
plot(EbNo_dB_List, blockErrorRate, 'LineWidth', 2);
title("Eb/N0 vs Block Error Rate (N,K) = (" + num2str(blockLength) + "," + num2str(msgLength) + ")");
xlabel('Eb/N0 (dB)');
ylabel('Block Error Rate');
grid on;