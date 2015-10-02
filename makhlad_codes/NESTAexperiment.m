N = 65536;
Sigma = 0.1;      %-- noise level
lambda = Sigma * sqrt(2 * log(N) );

Afor2f = @(signal) Hadamard2D_01(signal, M, 65536, Index_random_full(N / 2 + 1 : end, :) );
Aback2f = @(signal) Hadamard2Dtranspose_01(signal, M, 65536, Index_random_full(N / 2 + 1 : end, :) ).';
La = my_normest(Afor2f, Aback2f, N)^2;
nCalls = counter();
mu = .1 * Sigma / La; %--- can be chosen to be small 

opts = [];
opts.maxintiter = 6;
opts.TOlVar = 1e-6;
opts.verbose = 50;
opts.maxiter = 500;
opts.U = [];
opts.Ut = [];
opts.stoptest = 1;

[x, niter, resid, err, optsOut] = NESTA_UP(Afor2f, Aback2f, Y.', lambda, La, mu, opts);
x_NESTA = reshape(x, [], sqrt(N) );