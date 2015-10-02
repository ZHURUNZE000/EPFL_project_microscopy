%% function [x] = Aback2f(y)

function [x] = Aback2f(y)

MM = 256*2;
N = 256*256;

load('Index_random_full.mat');

signal = y;
x = Hadamard2Dtranspose_01(signal,MM,N,Index_random_full);
x = x';

end