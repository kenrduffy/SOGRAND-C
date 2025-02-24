% Guessing Random Additive Noise Decoding (GRAND)
% All code is subject to license:
% GRAND Codebase Non-Commercial Academic Research Use License 021722.pdf

% Simulation assesses the use of of blockwise SO from SOGRAND to manage
% undetected error rate.

% 1-line ORBGRAND 
% K. R. Duffy, W. An, and M. Medard, "Ordered reliability bits guessing 
% random additive noise decoding," IEEE Transactions on Signal Processing, 
% 70, 4528–4542, 2022.
% SISO SOGRAND  
% K. Galligan, P. Yuan, M. Médard & K. R. Duffy. Upgrade error detection 
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

clear;
%%
% Code parameters
n           = 31;
k           = 25;
code_class  = 'CRC';
[G, H]      = getGH_sys_CRC(n, k); % CRC from Koopman's database

%% Alternative classic code suggestion
% n           = 32;
% k           = 26;
% code_class  = 'eBCH';
% [G, H]      = getGH_sys_eBCH(n, k); % Extended BCH

%% Decoder parameters
L           = 2;
Tmax        = Inf;
thres       = Inf;
p_e         = 0.05;
%% Monte-Carlo parameters
EbN0dB      = 1:0.5:6;
NoErrors    = 20/p_e;
maxIt       = 10^9;
minIt       = 10^4;

%% Check if even code
if isequal(mod(sum(G,2),2),zeros(k,1))
    even = 1;
else
    even = 0;
end
%% Code and channel
EsN0dB      = EbN0dB + 10 * log10(2*k/n);
numSNR      = length(EsN0dB);

%% Loop over SNRs
BLER        = ones(1, numSNR);
UER         = ones(1, numSNR);
ER          = ones(1, numSNR);
MDR         = ones(1, numSNR);
NGavg       = zeros(1, numSNR);
for sp = 1:numSNR
    BlockError              = 0;
    undetectedBlockError    = 0;
    Erasure                 = 0;
    NG                      = 0;
    ntx                     = 0;
    scal = sqrt(10^(EsN0dB(sp) / 10));

    while ((BlockError < NoErrors && ntx < maxIt) || ntx < minIt)
        ntx = ntx + 1;
        msg = randsrc(k, 1, [0, 1]);
        c = mod(msg'*G, 2)';
        x = (1 - 2 * c) * scal;
        y = x + randn([n, 1]);
        llr = 2 * scal * y;

        [chat, N_guess, p_correct] = SOGRAND_blkSO(llr, H, L, Tmax, thres, even);
        NG = NG + N_guess;

        if p_correct > p_e
            flag = 0;
        else
            flag = 1;
        end

        if flag == 0
            BlockError = BlockError + 1;
            Erasure = Erasure + 1;
        elseif (~isequal(c, chat))
            BlockError = BlockError + 1;
            undetectedBlockError = undetectedBlockError + 1;
        end
    end
    BLER(sp)    = BlockError / ntx;
    UER(sp)   = undetectedBlockError / ntx;
    ER(sp)      = Erasure / ntx;
    MDR(sp)     = undetectedBlockError / BlockError;
    NGavg(sp)   = NG / (ntx);
    disp(['---' code_class '[' num2str(n) ',' num2str(k) ']---'])
    disp([num2str(EbN0dB(sp)), 'dB: BLER             = ', num2str(BLER(sp))]);
    disp([num2str(EbN0dB(sp)), 'dB: UER              = ', num2str(UER(sp))]);
    disp([num2str(EbN0dB(sp)), 'dB: MDR              = ', num2str(MDR(sp))]);
    disp([num2str(EbN0dB(sp)), 'dB: NGavg            = ', num2str(NGavg(sp))]);
    disp([num2str(EbN0dB(sp)), 'dB: NGavg/(info bit) = ', num2str(NGavg(sp)/k)]);
    save(['./results/uer-sogrand-' code_class '-' num2str(n) '-' num2str(k) '-' num2str(L) '-' num2str(p_e) '.mat'], 'EbN0dB', 'BLER', 'UER', 'NGavg');

end
