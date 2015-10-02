%% function [y] = Afor2f(x)

function [y] = Afor2f(x)

MM = 256*2;
N = 256*256;

load('Index_random_full.mat');

x = reshape(x,[],1);
signal = x;
y = Hadamard2D_01(signal,MM,N,Index_random_full);

end