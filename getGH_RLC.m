function [G,H] = getGH_RLC(n, k)
P = randi([0 1], k, n-k);
G = [eye(k), P];
H = [P', eye(n-k)];
end