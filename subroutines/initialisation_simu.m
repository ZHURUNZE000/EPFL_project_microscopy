%% prior dependent quantities
param_1 = 0; 
param_2 = 0;
param_3 = 0;
switch opt.prior    
    case 'L1'
        param_1 = opt.L1_min; 
        param_2 = opt.L1_max;
    case 'GaussExponential'
        param_1 = opt.GaussExponential_expo; 
        param_2 = opt.GaussExponential_mGauss; 
        param_3 = opt.GaussExponential_varGauss;
    case 'GaussLaplace'
        param_1 = opt.GaussLaplace_expo; 
        param_2 = opt.GaussLaplace_mGauss; 
        param_3 = opt.GaussLaplace_varGauss;
    case 'SparseExponential'
        param_1 = opt.SparseExponential_expo;
end

%% size stuff
M = opt.M; 
N = opt.N;
alpha_ = M / N;

%% remove mean stuff
Gmean = HT(ones(1, M) ) / M; 	
% Gmean = min(Gmean, 0.9);
Gm2 = Gmean.^2;
Yeff = Y - mean(Y);

%% some inputs
if opt.signal_rho > 0
    rho = opt.signal_rho;   
else
    rho = alpha_ / 10;
end
varNoise = opt.varNoise;

%% AMP fields initialisation 
W = Yeff;
V = ones(1, M);
R = zeros(1, N); 
S2 = zeros(1, N);
av_mess = zeros(1, N);
var_mess = rho * ones(1, N);