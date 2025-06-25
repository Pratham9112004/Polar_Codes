% Polar Code Encoding and Successive Cancellation List Decoding 

clearvars;

load('Reliability_Sequence.mat');
blockLength = 512;
n = log2(blockLength);
A = 249;
crcBits = 11;
crcPoly = fliplr([1 1 1 0 0 0 1 0 0 0 0 1]);
infoBits = A + crcBits;
codeRate = A / blockLength;

global L_size;
L_size = 4;

reliableSet = Reliability_Sequence(Reliability_Sequence <= blockLength);
frozenSet = reliableSet(1:blockLength - infoBits);

NSim = 400;
EbNo_dB_List = 0:0.5:10;

successRate = zeros(1, length(EbNo_dB_List));
bitErrorRate = zeros(1, length(EbNo_dB_List));
blockErrorRate = zeros(1, length(EbNo_dB_List));

for idx = 1:length(EbNo_dB_List)
    Eb_dB = EbNo_dB_List(idx);
    successfulBlocks = 0;
    bitErrorsTotal = 0;
    blockErrorsTotal = 0;

    for trial = 1:NSim
        msg_input = randi([0 1], 1, A);
        [~, rem] = gfdeconv([zeros(1, crcBits) fliplr(msg_input)], crcPoly);
        paddedData = [msg_input fliplr([rem zeros(1, crcBits - length(rem))])];
        msgVec = zeros(1, blockLength);
        msgVec(reliableSet(blockLength - infoBits + 1:end)) = paddedData;

        encodedVec = polarEncode(msgVec, n);
        rxVec = addBPSKandAWGN(encodedVec, Eb_dB, A, blockLength);

        global LogL;
        global decisions;
        global pathMetrics;
        global nodeTracker;

        LogL = zeros(L_size, n+1, blockLength);
        LogL(:,1,:) = repmat(rxVec, L_size, 1, 1);
        decisions = zeros(L_size, n+1, blockLength);
        nodeTracker = zeros(1, 2 * blockLength - 1);
        pathMetrics = Inf * ones(L_size, 1);
        pathMetrics(1) = 0;

        listDecoder(0, 0, n, frozenSet);

        decodedSet = squeeze(decisions(:, n+1, reliableSet(blockLength - infoBits + 1:end)));
        chosen = 1;
        for l = 1:L_size
            [~, crcR] = gfdeconv(fliplr(decodedSet(l,:)), crcPoly);
            if isequal(crcR, 0)
                chosen = l;
                break;
            end
        end

        message_cap = decodedSet(chosen, 1:A);

        currErrors = sum(msg_input ~= message_cap);
        if currErrors == 0
            successfulBlocks = successfulBlocks + 1;
        else
            bitErrorsTotal = bitErrorsTotal + currErrors;
            blockErrorsTotal = blockErrorsTotal + 1;
        end
    end

    successRate(idx) = successfulBlocks / NSim;
    bitErrorRate(idx) = bitErrorsTotal / (A * NSim);
    blockErrorRate(idx) = blockErrorsTotal / NSim;
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
    noiseStd = sqrt(1 / (2 * Es));
    noise = noiseStd * randn(1, N);
    y = bpsk + noise;
end

function listDecoder(depth, node, n, frozenSet)
    global LogL;
    global decisions;
    global pathMetrics;
    global L_size;
    global nodeTracker;

    fOp = @(a,b) sign(a).*sign(b).*min(abs(a), abs(b));
    gOp = @(a,b,c) b + (1 - 2 * c) .* a;

    pos = 2^depth - 1 + node + 1;

    if depth == n
        nodeTracker(pos) = 3;
        Lvals = squeeze(LogL(:, depth+1, node+1));
        if any(frozenSet == node+1)
            decisions(:, depth+1, node+1) = 0;
            pathMetrics = pathMetrics + abs(Lvals).*(Lvals < 0);
        else
            bitDec = Lvals < 0;
            newPaths = [pathMetrics; pathMetrics + abs(Lvals)];
            [pathMetrics, idx] = mink(newPaths, L_size);
            keep = idx > L_size;
            idx(keep) = idx(keep) - L_size;
            bitDec = bitDec(idx);
            bitDec(keep) = 1 - bitDec(keep);
            LogL = LogL(idx, :, :);
            decisions = decisions(idx, :, :);
            decisions(:, depth+1, node+1) = bitDec;
        end
        return;
    end

    segLen = 2^(n - depth);
    r = squeeze(LogL(:, depth+1, segLen*node+1 : segLen*(node+1)));

    if nodeTracker(pos) == 0
        nodeTracker(pos) = 1;
        a = r(:, 1:segLen/2);
        b = r(:, segLen/2+1:end);
        LogL(:, depth+2, segLen*2*node+1 : segLen*(2*node+1)) = fOp(a,b);
        listDecoder(depth+1, 2*node, n, frozenSet);
    end

    r = squeeze(LogL(:, depth+1, segLen*node+1 : segLen*(node+1)));
    if nodeTracker(pos) == 1
        nodeTracker(pos) = 2;
        a = r(:, 1:segLen/2);
        b = r(:, segLen/2+1:end);
        uL = squeeze(decisions(:, depth+2, segLen*2*node+1 : segLen*(2*node+1)));
        LogL(:, depth+2, segLen*(2*node+1)+1 : segLen*(2*node+2)) = gOp(a, b, uL);
        listDecoder(depth+1, 2*node+1, n, frozenSet);
    end

    if nodeTracker(pos) == 2
        nodeTracker(pos) = 3;
        uL = squeeze(decisions(:, depth+2, segLen*2*node+1 : segLen*(2*node+1)));
        uR = squeeze(decisions(:, depth+2, segLen*(2*node+1)+1 : segLen*(2*node+2)));
        decisions(:, depth+1, segLen*node+1 : segLen*(node+1)) = ...
            [mod(uL + uR, 2), uR];
    end
end

figure;
plot(EbNo_dB_List, successRate, 'LineWidth', 2);
title("Success Rate vs Eb/N0  (N,K)=(" + blockLength + "," + infoBits + ")");
xlabel('Eb/N0 (dB)');
ylabel('Success Rate');
grid on;

figure;
plot(EbNo_dB_List, bitErrorRate, 'LineWidth', 2);
title("Bit Error Rate vs Eb/N0  (N,K)=(" + blockLength + "," + infoBits + ")");
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Rate');
grid on;

figure;
plot(EbNo_dB_List, blockErrorRate, 'LineWidth', 2);
title("Block Error Rate vs Eb/N0  (N,K)=(" + blockLength + "," + infoBits + ")");
xlabel('Eb/N0 (dB)');
ylabel('Block Error Rate');
grid on;