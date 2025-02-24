% Guessing Random Additive Noise Decoding (GRAND)
% All code is subject to license:
% GRAND Codebase Non-Commercial Academic Research Use License 021722.pdf

function [L_APP, L_E, N_guess] = SOGRAND_bitSO(llr, H, L, Tmax, thres, even)
    [chat_list, s_list, N_guess, curL, pNL, ~] = SOGRAND_mex(llr, uint8(reshape(H', [], 1)), int32(-1), uint64(L), uint64(Tmax), thres, even);
    
    chat = chat_list(:,1:curL);
    PM = s_list(2, 1:curL);
    p_notinlist = max(pNL, eps);
    
    pp1 = 1./(1+exp(llr));
    pp0 = 1 - pp1;
    
    pp1 = max(pp1, eps); pp1 = min(pp1, 1-eps);
    pp0 = max(pp0, eps); pp0 = min(pp0, 1-eps);
    
    p = exp(-PM);
    
    p1 = sum(chat.*p,2);
    p0 = sum((1-chat).*p,2);
    
    p0 = p0 + p_notinlist * pp0;
    p1 = p1 + p_notinlist * pp1;
    
    L_APP = log(p0) - log(p1);
    L_E = L_APP - llr;
end