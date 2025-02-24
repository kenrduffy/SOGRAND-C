% Guessing Random Additive Noise Decoding (GRAND)
% All code is subject to license:
% GRAND Codebase Non-Commercial Academic Research Use License 021722.pdf

% Simulation assesses the calibration of blockwise SO from SOGRAND.

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
%% Monte-Carlo parameters
EbN0dB = 2;
maxIt = 10^5;
minlogp = 5;
step = 0.5;
r = 10.^-(0:step:minlogp);
%% Code parameters
n = 32;
k = 26;
code_class = 'RLC';
[G, H] = getGH_RLC(n, k); % Random linear code

%% Alternative CRC code suggestion
% n           = 31;
% k           = 25;
% code_class  = 'CRC';
% [G, H]      = getGH_sys_CRC(n, k); % CRC from Koopman's database

%% Alternative classic TPC code suggestion
% n           = 32;
% k           = 26;
% code_class  = 'eBCH';
% [G, H]      = getGH_sys_eBCH(n, k); % Extended BCH

%% Check even code
if isequal(mod(sum(G,2),2),zeros(k,1))
    even = 1;
else
    even = 0;
end

%% Decoder parameters
L = 4;
Tmax = Inf;
%% 
R = k / n;
EsN0dB = EbN0dB + 10 * log10(2*k/n);
scal = sqrt(10^(EsN0dB / 10));
%% Loop over SNRs
n_sample = zeros(1, length(r));
n_negative = zeros(1, length(r));
n_p = zeros(1, length(r));
for ntx = 1:maxIt
    msg = randsrc(k, 1, [0, 1]);
    c = mod(msg'*G, 2)';
    %c = zeros(n, 1);

    x = (1 - 2 * c) * scal;
    y = x + randn([n, 1]);
    llr = 2 * scal * y;

    % confidence
    [chat, N_guess, p_correct, chat_list, p_inList] = SOGRAND_blkSO(llr, H, L, Tmax, Inf, even);

    ind = ceil(-(log10(p_inList) / step));
    ind = max(1, ind);
    ind = min(ind, length(r));
    n_p(ind) = n_p(ind) + p_inList;
    n_sample(ind) = n_sample(ind) + 1;

    %
    flag = 0;
    for l = 1:L
        if isequal(chat_list(:, l), c)
            flag = 1;
            break
        end
    end
    %
    if flag == 0
        n_negative(ind) = n_negative(ind) + 1;
    end
    %
    if mod(ntx, maxIt/10) == 0
        disp(['@ ', num2str(ntx), ' / ', num2str(maxIt)])
    end
end

figure(1)
clf
loglog([10^(-minlogp - 1), 1], [10^(-minlogp - 1), 1], 'k','LineWidth',2)
hold on
grid on
p = n_negative ./ n_sample;
pp = n_p ./ n_sample;
loglog(pp, p, 'rx-','LineWidth',2)
grid 'on'
xlabel('Predicted')
ylabel('Conditional')
set(gca,'FontSize',16)
