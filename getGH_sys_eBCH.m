function [G, H] = getGH_sys_eBCH(n, k)
    [~, H] = getGH_BCH(n-1, k);
    
    H = [H, zeros(n-k-1, 1)];
    H = [H; ones(1, n)];
    for i = 1:n - k - 1
        H(n-k, :) = mod(H(i, :)+H(n-k, :), 2);
    end
    
    P = H(:, 1:k);
    G = [eye(k), P'];
end

function [G, H] = getGH_BCH(n, k)
    P = [];
    for i = 1:k
        w = zeros(k, 1);
        w(i) = 1;
        p = bchenc(gf(w'), n, k);
        P = [P; p.x(k+1:n)];
    end
    P = double(P);
    G = [eye(k), P];
    H = [P', eye(n-k)];
end
