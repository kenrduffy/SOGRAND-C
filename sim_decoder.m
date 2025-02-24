% Guessing Random Additive Noise Decoding (GRAND)
% All code is subject to license:
% GRAND Codebase Non-Commercial Academic Research Use License 021722.pdf

% Simulation peforms decoding with basic or 1-line OBRGRAND.

% 1-line ORBGRAND 
% K. R. Duffy, W. An, and M. Medard, "Ordered reliability bits guessing 
% random additive noise decoding," IEEE Transactions on Signal Processing, 
% 70, 4528–4542, 2022.

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

clear;
%% Monte-Carlo parameters
EbN0dB      = 4:0.5:6;
NoErrors    = 50;
maxIt       = 10^6;
minIt       = 10^2;
% Code parameters
n           = 129;
k           = 114;
code_class  = 'CRC';
[G, H]      = getGH_sys_CRC(n, k); % CRC from Koopman's database

%% Alternative classic code suggestion
% n           = 128;
% k           = 113;
% code_class  = 'eBCH';
% [G, H]      = getGH_sys_eBCH(n, k); % Extended BCH

%% Decoder parameters
IC          = -1;  % Set to 0 for basic ORBGRAND
L           = 1;   % Maximum list size
Tmax        = Inf; % Maximum number of queries per
thres       = 1;   % Abandon decoding if list already has >thres prob. 

%% Check even code
if isequal(mod(sum(G, 2), 2), zeros(k, 1))
    even = 1;
else
    even = 0;
end
%% Code and channel
R           = (k / n);
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
    sigma = 1 / sqrt(10^(EsN0dB(sp) / 10));
    while ((BlockError < NoErrors && ntx < maxIt) || ntx < minIt)
        ntx = ntx + 1;
        c = zeros(1, n); % To store codeword.
        %% Enc
        u = randsrc(1, k, [0, 1]);
        c= mod(u*G, 2); 
        %% binary input AWGN channel
        x = 1 - 2 * c;
        y = x + sigma * randn([1, n]);
        L_channel = 2 * y / (sigma^2);
        %% Decoding
        [c_HD, ~, N_guess, ~, ~, ~] = SOGRAND_mex(L_channel', uint8(reshape(H', [], 1)), int32(IC), uint64(L), uint64(Tmax), 0, even);
        c_HD =c_HD';
        NG = NG + N_guess;
        %% error collection
        if (~isequal(c, c_HD))
            BlockError = BlockError + 1;
            uhat = c_HD(1, 1:k);
            BitError = BitError + sum(uhat(:) ~= u(:));
        end     
    end
    disp(['---' code_class '[' num2str(n) ',' num2str(k) ']---Eb/N0 dB ', num2str(EbN0dB(sp)), ' dB:---'])
    BLER(sp)    = BlockError / ntx;
    BER(sp)     = BitError / (ntx * k);
    Iavg(sp)    = n_iter / ntx;
    NGavg(sp)   = NG / (ntx);
    disp([' BLER             = ', num2str(BLER(sp))]);
    disp([' BER              = ', num2str(BER(sp))]);
    disp([' NGavg            = ', num2str(NGavg(sp))]);    
    disp([' NGavg/(info bit) = ', num2str(NGavg(sp)/k)]);

    save(['./results/sogrand-' code_class '-' num2str(n) '-' num2str(k) '.mat'], 'EbN0dB', 'BLER', 'BER', 'Iavg', 'NGavg', 'NGavg_p', 'thres', 'G', 'H');
end
