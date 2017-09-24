function [ trans] = GFT(S, x, k)
%GFT Summary of this function goes here
%   Detailed explanation goes here

[V, ~] = eigs(S, k);

H = V';
trans = H*x;

end

