% Guessing Random Additive Noise Decoding (GRAND)
% All code is subject to license:
% GRAND Codebase Non-Commercial Academic Research Use License 021722.pdf

function [chat, N_guess, p_correct, chat_list, p_inList] = SOGRAND_blkSO(llr, H, L, Tmax, thres, even)
    [chat_list, s_list, N_guess, curL, pNL, ~] = SOGRAND_mex(llr, uint8(reshape(H', [], 1)), int32(-1), uint64(L), uint64(Tmax), thres, even);
    pNL = max(eps, pNL);
    i_ml = find(s_list(2, 1:curL) == min(s_list(2, 1:curL)));
    chat = chat_list(:, i_ml);
    pC = pNL + sum(exp(-s_list(2, :)));
    p_correct = 1 - (exp(-s_list(2, i_ml)) / pC);
    p_inList = 1 - s_list(4,curL);
end
