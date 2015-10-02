close all; clear; clc;

load Index_random_full.mat

%% number of measurements
alpha_ = 0.01;
N = 256^2;
M = round(alpha_ * N);

%% measure and image
sparsity = 5e-5;
noise_ = 0;
sizeBeads = 2;
intensity = 2000;
% F = (randn(M, N) > 0);
% Afor2f = @(x) F * x.';
% Aback2f = @(x) x * F;
% Aback2f = @(x) Hadamard2Dtranspose_01(x, M, N, Index_random_full); 
% Afor2f = @(x) Hadamard2D_01(x, M, N, Index_random_full);



numBlockC = 4;
numBlockL = 4;
JJ = .1;
density = .1;
w = 1;
alphaCs = zeros(1, numBlockL); 
alphaCs(1) = alpha_ * 1;
if (numBlockC > 1); alphaCs(2 : numBlockL) = (numBlockC .* alpha_ - alphaCs(1) ) ./ (numBlockL - 1); end % measurement rate of the bulk blocks, taken into acount if numBlockC > 1 : not to be modified
Nblock = N / numBlockC; % number of variables per block in the seeded matrix, must divide N (POWER OF 2 if type = 'real')
for i = 1 : numBlockL; Mblock(i) = floor(alphaCs(i) * Nblock); end % number of lines of the blocks in the seeded matrix
Mblock = min(Mblock, Nblock - 1);
M = sum(Mblock);
J = createSeededJ(numBlockL, numBlockC, JJ, w, 1, 0);
F = createSeededRandomMatrix_density(J, density, Mblock, Nblock);
Afor2f = @(x) F * x(:);
Aback2f = @(x) x(:)' * F;



[original, Y, Y0, noise] = CFM_sim_simple(Afor2f, sparsity, intensity, sizeBeads, noise_);
% resc = sqrt(mean(Y.^2) );
% Y = Y' / resc;
% Y0 = Y0' / resc;
% F = F / resc;

Y = Y';

for i = 1 : M
	el = std(F(i, :) );
	F(i , :) = F(i , :) ./ el;
	Y(i) = Y(i) ./ el;
end

% Aback2f = @(x) Hadamard2Dtranspose_01(x, M, N, Index_random_full) / resc; 
% Afor2f = @(x) Hadamard2D_01(x, M, N, Index_random_full) / resc;
% F = F / N;
Afor2f = @(x) F * x.';
Aback2f = @(x) x * F;

%% AMP parameters
postProcess = 1; % Do you want to post process to show only the pixels that are not noise from the AMP point of view ?
opt.prior = 'GaussExponential'; % 'L1', 'GaussExponential', 'SparseExponential' or 'GaussLaplace'
opt.tMax = 500;
opt.print = 1;
opt.conv_ = 1e-7;
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
opt.N = N;
opt.signal_rho = sparsity * 10;
opt.L1_min = -Inf; 
opt.L1_max = Inf;
opt.GaussExponential_expo = 3; 
opt.GaussExponential_mGauss = 0; 
opt.GaussExponential_varGauss = 1e-3;
opt.SparseExponential_expo = 3;
opt.GaussLaplace_expo = 1; 
opt.GaussLaplace_mGauss = 0; 
opt.GaussLaplace_varGauss = 0.001;
opt.showImage = 1;
opt.weightMf = 0;

%% AMP reconstruction
disp('AMP reconstruction');
if strcmp(opt.prior, 'GaussExponential') || strcmp(opt.prior, 'SparseExponential') || strcmp(opt.prior, 'GaussLaplace')
	[X, weightNoise] = CSBP_Solver_simu(Y, Afor2f, Aback2f, opt);
else
	X = CSBP_Solver_simu(Y, Afor2f, Aback2f, opt);
end

% %% FastIHT
% disp('FastIHT reconstruction');
% x_FastIHT = Fast_IHT_v2(256, Y.', Afor2f, Aback2f, 1, 500, 3, 1);
% x_FastIHT = x_FastIHT(:);

% %% NESTA;
% disp('NESTA reconstruction');
% NESTAexperiment;
% x_NESTA = x_NESTA(:);

%% post processing amplification procedure for AMP
if postProcess && (strcmp(opt.prior, 'GaussExponential') || strcmp(opt.prior, 'SparseExponential') || strcmp(opt.prior, 'GaussLaplace') )
	X(weightNoise > 0.5) = 0;
	% X = amplifyIntensities(X, 1);
	ampli = fspecial('disk', 2); % filter to make the beads as disks
	X = imfilter(reshape(X, 256, [] ), ampli);

	% trueRho = mean(original(:) > 8000);	
	% numStd = norminv(1 - trueRho / 2, 0, 1);
	 
	% x_NESTA(x_NESTA < mean(x_NESTA) + numStd * std(x_NESTA(:) ) ) = 0;
	% x_NESTA = amplifyIntensities(x_NESTA, 1);

	% x_FastIHT(x_FastIHT < mean(x_FastIHT) + 1 * std(x_FastIHT) ) = 0;
	% x_FastIHT = amplifyIntensities(x_FastIHT, 1);	
end

% show images
subplot(2, 2, 1); imshow(reshape(original, 256, [] ) ); title('original');
subplot(2, 2, 2); imshow(reshape(X, 256, [] ) ); title(['AMP M=', num2str(M) ] );
% subplot(2, 2, 3); imagesc(reshape(x_FastIHT, 256, [] ) ); title(['Fast IHT M=', num2str(M) ] );
% subplot(2, 2, 4); imagesc(reshape(x_NESTA, 256, [] ) ); title(['NESTA M=', num2str(M) ] );