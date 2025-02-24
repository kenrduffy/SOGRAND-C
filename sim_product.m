% Guessing Random Additive Noise Decoding (GRAND)
% All code is subject to license:
% GRAND Codebase Non-Commercial Academic Research Use License 021722.pdf

% Simulation peforms turbo product code decoding with SOGRAND as the SISO
% component decoder.

% 1-line ORBGRAND 
% K. R. Duffy, W. An, and M. Medard, "Ordered reliability bits guessing 
% random additive noise decoding," IEEE Transactions on Signal Processing, 
% 70, 4528–4542, 2022.
% SISO SOGRAND  
% K. Galligan, P. Yuan, M. Médard & K. R. Duffy. "Upgrade error detection 
% to prediction with GRAND". Proceedings of Globecom, 2023.
% P. Yuan, M. Médard, K. Galligan & K. R. Duffy, "Soft-output (SO) GRAND 
% and long, low rate codes to outperform 5 LDPC codes", IEEE Transactions
% on Wireless Communications, 2025.

% Exploiting codebook structure to alter query order in GRAND algorithms 
% was first proposed in 
% M. Rowshan and J. Yuan, "Constrained Error Pattern Generation for GRAND", 
% IEEE International Symposium on Information Theory, 2022, and expanded in
% M. Rowshan and J. Yuan, "Segmented GRAND: Combining sub-patterns in 
% near-ML order,” IEEE Transactions on Communications, 2025.
% In particular, for even codes, ORBGRAND algorithms implemented with 
% landslide algorithms can skip queries without creating them with no change
% in decoding performance. That is exploited in ORBGRAND here.

% The C implementation of the component decoder here DOES NOT avail of the 
% parallelizability of querying and nor does it pipeline sorting 
% reliabilities with querying. ASIC implementations avail of those 
% features to reduce energy and latency.
% A. Riaz, A. Yasar, F. Ercan, W. An, J. Ngo, K. Galligan, M. Médard, 
% K. R. Duffy & R. T. Yazicigil. "A Sub-0.8-pJ/bit Universal Soft-Detection 
% Decoder Using ORBGRAND." IEEE Journal of Solid State Circuits, 2025.

% Product code decoding is also highly parallelizable, where all rows can 
% be decoded in parallel and all columns can be decoded in parallel. That 
% is not exploited in this implementation, which is purely for
% illustration.

% LLR updates are executed in MATLAB rather than C.

clear;
%% Monte-Carlo parameters
EbN0dB      = 1.5:0.25:2.5;
NoErrors    = 20;
maxIt       = 10^6;
minIt       = 10^2;
%% Code parameters
n           = 31;
k           = 25;
code_class  = 'CRC';
[G, H]      = getGH_sys_CRC(n, k); % CRC from Koopman's database

%% Alternative classic TPC component code suggestion
% n           = 32;
% k           = 26;
% code_class  = 'eBCH';
% [G, H]      = getGH_sys_eBCH(n, k); % Extended BCH

%% Decoder parameters
L           = 4;  % Maximum list size
Imax        = 20; % maximum number of iterations
Tmax        = Inf; % Maximum number of queries per componet decoding
p_ET        = 1e-5; 
thres       = 1 - p_ET; % Abandon decoding if list already has >thres prob. 
alpha       = 0.5 * ones(1, 50); % Extrinsic LLR scaling

%% Check even code
if isequal(mod(sum(G, 2), 2), zeros(k, 1))
    even = 1;
else
    even = 0;
end
%% Code and channel
R           = (k / n)^2;
EsN0dB      = EbN0dB + 10 * log10(2*R);
numSNR      = length(EsN0dB);

%% Loop over SNRs
BLER        = zeros(1, numSNR);
BER         = zeros(1, numSNR);
Iavg        = zeros(1, numSNR);
NGavg       = zeros(1, numSNR);
NGavg_p     = zeros(1, numSNR);

for sp = 1:numSNR
    BlockError  = 0;
    BitError    = 0;
    n_iter      = 0;
    ntx         = 0;
    NG          = 0;
    NG_p        = 0;
    sigma = 1 / sqrt(10^(EsN0dB(sp) / 10));
    while ((BlockError < NoErrors && ntx < maxIt) || ntx < minIt)
        ntx = ntx + 1;
        c = zeros(n, n); % To store codeword.
        %% Enc
        u = randsrc(k, k, [0, 1]);
        for row = 1:k
            c(row, :) = mod(u(row, :)*G, 2);
        end
        for col = 1:n
            c(:, col) = mod(c(1:k, col)'*G, 2);
        end
        %% binary input AWGN channel
        x = 1 - 2 * c;
        y = x + sigma * randn([n, n]);
        L_channel = 2 * y / (sigma^2);
        %% Decoding
        L_APP = zeros(size(L_channel));
        L_E = zeros(size(L_channel));
        L_A = zeros(size(L_channel));
        c_HD = 0.5 * (1 - sign(L_channel));
        for i = 1:Imax
            %% row
            NGmax = 0;
            n_iter = n_iter + 0.5;
            L_A = alpha(2*i-1) * L_E;
            input = L_channel + L_A;
            for row = 1:n
                [L_APP(row, :), L_E(row, :), N_guess] = SOGRAND_bitSO(input(row, :)', H, L, Tmax, thres, even);
                NG = NG + N_guess;
                NGmax = max(NGmax, N_guess);
            end
            NG_p = NG_p + NGmax;
            % early-termination
            c_HD = 0.5 * (1 - sign(L_APP));
            s1 = mod(c_HD*H', 2);
            s2 = mod(c_HD'*H', 2);
            if sum(s1(:)) == 0 && sum(s2(:)) == 0
                break
            end
            %% column
            NGmax = 0;
            n_iter = n_iter + 0.5;
            L_A = alpha(2*i) * L_E;
            input = L_channel + L_A;
            for col = 1:n
                [L_APP(:, col), L_E(:, col), N_guess] = SOGRAND_bitSO(input(:, col), H, L, Tmax, thres, even);
                NG = NG + N_guess;
                NGmax = max(NGmax, N_guess);
            end
            NG_p = NG_p + NGmax;
            %% early-termination
            c_HD = 0.5 * (1 - sign(L_APP));
            s1 = mod(c_HD*H', 2);
            s2 = mod(c_HD'*H', 2);
            if sum(s1(:)) == 0 && sum(s2(:)) == 0
                break
            end
        end
        %% error collection
        if (~isequal(c, c_HD))
            BlockError = BlockError + 1;
            uhat = c_HD(1:k, 1:k);
            BitError = BitError + sum(uhat(:) ~= u(:));
        end
        
    end
    disp(['---' code_class '^2, [' num2str(n^2) ',' num2str(k^2) ']---Eb/N0 dB ', num2str(EbN0dB(sp)), ' dB:---'])
    BLER(sp)    = BlockError / ntx;
    BER(sp)     = BitError / (ntx * k * k);
    Iavg(sp)    = n_iter / ntx;
    NGavg(sp)   = NG / (ntx);
    NGavg_p(sp) = NG_p / (ntx);
    disp([' BLER             = ', num2str(BLER(sp))]);
    disp([' BER              = ', num2str(BER(sp))]);
    disp([' Iavg             = ', num2str(Iavg(sp))]);
    disp([' NGavg            = ', num2str(NGavg(sp))]);    
    disp([' NGavg/(info bit) = ', num2str(NGavg(sp)/k^2)]);
    disp([' NGavg_p          = ', num2str(NGavg_p(sp))]);

    save(['./results/prod-sogrand-' code_class '-' num2str(n) '-' num2str(k) '-' num2str(L) '-' num2str(-log10(p_ET)) '.mat'], 'EbN0dB', 'BLER', 'BER', 'Iavg', 'NGavg', 'NGavg_p', 'thres', 'G', 'H');
end
