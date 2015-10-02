%% data
load Index_random_full.mat
load acquisition_part1.mat
original = imread('snap_part1.tif'); % the image to compare with, i.e. the true one

%% number of measurements
M = 512;

%% measure
Y = data2b(1 : M);

%% the operators to use WARNING: non zero mean operators (elements are 0,1)
H = @(x) Hadamard2D_01(x, M, 65536, Index_random_full); % forward op
HT = @(x) Hadamard2Dtranspose_01(x, M, 65536, Index_random_full); % backward op

%% put you solver here :)
% results = yourSolver(Y, H, HT);