%% WORKS ONLY FOR BINARY 0/1 MATRICES

%% AMP computation of the gaussian fields
V_new = H(prior.var_mess).' + Gm2 * prior.var_mess.' - 2 * H(prior.var_mess .* Gmean).';
W_new = dumping(W, (H(prior.av_mess) - Gmean * prior.av_mess.').' - ((Yeff - W) ./ (varNoise + V) ) .* V_new, opt.dump_mess);
V_new = dumping(V, V_new, opt.dump_mess);
PRELI = (varNoise + V_new).^(-1);
PRELI2 = (Yeff - W_new) .* PRELI;
var_1 = max(1e-20, sum(PRELI) * Gm2 + HT(PRELI) .* (1 - 2 * Gmean) );
var_2 = HT(PRELI2) - sum(PRELI2) * Gmean;

% V_new = (G.^2 * prior.var_mess.' + Gm1d2 * prior.var_mess.' - (2 .* G .* (ones(M,1) * sum(G) ./ M) ) * prior.var_mess.').';
% W_new = dumping(W,(G * prior.av_mess.' - Gmean1d * prior.av_mess.').' - ((Yeff - W) ./ (n_and_e.var_noise + V) ) .* V_new,opt.dump_mes);
% V_new = dumping(V,V_new,opt.dump_mes);
% PRELI = 1 ./ (n_and_e.var_noise + V_new); % [1,M]
% PRELI2 = ((Yeff - W_new) ./ (n_and_e.var_noise + V_new) ); % [M,1]
% var_1 = max(1e-18,PRELI * G.^2 + sum(PRELI) * Gm1d2 - 2 .* PRELI * (G .* (ones(M,1) * sum(G) ./ M) ) );
% var_2 = PRELI2 * G - sum(PRELI2) * Gmean1d;

%% keep track of the previous step
av_mess_old = prior.av_mess;
var_mess_old = prior.var_mess;

%% update
if t > 1 && (strcmp(opt.prior, 'GaussExponential') || strcmp(opt.prior, 'SparseExponential') || strcmp(opt.prior, 'GaussLaplace') ) && (opt.weightMf > 0)
	computeMeanField;
	prior.R = var_2 ./ var_1 + prior.av_mess + mf * opt.weightMf;
else 
	prior.R = var_2 ./ var_1 + prior.av_mess; 	
end
prior.S2 = var_1.^(-1);

%% denoising with prior
prior = F(prior);
V = V_new;
W = W_new;