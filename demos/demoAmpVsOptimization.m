close all; clear; clc;

othersAlgo = 0;

%% data
load Index_random_full.mat
load acquisition_part1.mat
original = imread('snap_part1.tif');

%% number of measurements
M = 512; % for M = 512, 1024, 2048, 4096 or 8192, results are already in memory

%% measure
Y = data2b(1 : M);

%% AMP parameters
postProcess = 1; % Do you want to post process to show only the pixels that are not noise from the AMP point of view ?
opt.prior = 'GaussExponential'; % 'L1', 'GaussExponential', 'SparseExponential' or 'GaussLaplace'
opt.tMax = 50;
opt.print = 1;
opt.conv_ = 5e-4;
opt.learn = 1;
opt.learnNoise = 0;
opt.varNoise = 0;
if strcmp(opt.prior, 'GaussExponential') || strcmp(opt.prior, 'SparseExponential') || strcmp(opt.prior, 'GaussLaplace')
	opt.dump_mess = 0.2;
elseif strcmp(opt.prior, 'L1')
	opt.dump_mess = 0.9;
end
opt.dump_learn = 0.9;
opt.M = M;
opt.N = 256^2;
opt.signal_rho = 1e-3;
opt.L1_min = -Inf; 
opt.L1_max = Inf;
opt.GaussExponential_expo = 3; 
opt.GaussExponential_mGauss = 0; 
opt.GaussExponential_varGauss = 0.001;
opt.SparseExponential_expo = 3;
opt.GaussLaplace_expo = 1; 
opt.GaussLaplace_mGauss = 0; 
opt.GaussLaplace_varGauss = 0.001;
opt.showImage = 1;
opt.weightMf = .1;
opt.part2 = 0;

%% AMP reconstruction
disp('AMP reconstruction');
if strcmp(opt.prior, 'GaussExponential') || strcmp(opt.prior, 'SparseExponential') || strcmp(opt.prior, 'GaussLaplace')
	[X, weightNoise] = CSBP_Solver(Y, opt);
else
	X = CSBP_Solver(Y, opt);
end

%% experiment with FastIHT and NESTA
if othersAlgo
	if (M == 512) || (M == 1024) || (M == 2048) || (M == 4096) || (M == 8192)	
		load ResultsOptimization.mat
		x_NESTA = eval(['x_NESTA_', num2str(M) ] );
		x_FastIHT = eval(['x_FastIHT_', num2str(M) ] );
		x_NESTA = x_NESTA(:);
		x_FastIHT = x_FastIHT(:);
	else
		%% FastIHT
		disp('FastIHT reconstruction');
		Afor2f = @(signal) Hadamard2D_01(signal, M, 65536, Index_random_full);
		Aback2f = @(signal) Hadamard2Dtranspose_01(signal, M, 65536, Index_random_full);
		x_FastIHT = Fast_IHT_v2(256, Y.', Afor2f, Aback2f, 1, 500, 3, 1);
		x_FastIHT = x_FastIHT(:);

		%% NESTA;
		disp('NESTA reconstruction');
		NESTAexperiment;
		x_NESTA = x_NESTA(:);
	end
end

%% post processing amplification procedure for AMP
if postProcess && (strcmp(opt.prior, 'GaussExponential') || strcmp(opt.prior, 'SparseExponential') || strcmp(opt.prior, 'GaussLaplace') )
	X(weightNoise > 0.5) = 0;
	X = amplifyIntensities(X, 1);

	if othersAlgo
		trueRho = mean(original(:) > 8000);	
		numStd = norminv(1 - trueRho / 2, 0, 1);
		 
		x_NESTA(x_NESTA < mean(x_NESTA) + numStd * std(x_NESTA(:) ) ) = 0;
		x_NESTA = amplifyIntensities(x_NESTA, 1);

		x_FastIHT(x_FastIHT < mean(x_FastIHT) + numStd * std(x_FastIHT) ) = 0;
		x_FastIHT = amplifyIntensities(x_FastIHT, 1);	
	end
end

% show images
subplot(2, 2, 1); imagesc(original); title('original');
subplot(2, 2, 2); imagesc(reshape(X, 256, [] ) ); title(['AMP M=', num2str(M) ] );
if othersAlgo
	subplot(2, 2, 3); imagesc(reshape(x_FastIHT, 256, [] ) ); title(['Fast IHT M=', num2str(M) ] );
	subplot(2, 2, 4); imagesc(reshape(x_NESTA, 256, [] ) ); title(['NESTA M=', num2str(M) ] );
end