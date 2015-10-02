%% AMP computation of the gaussian fields
V_new = H(prior.var_mess).' + Gm2 * prior.var_mess.' - 2 * H(prior.var_mess .* Gmean).';
W_new = dumping(W, (H(prior.av_mess) - Gmean * prior.av_mess.').' - ((Yeff - W) ./ (varNoise + V) ) .* V_new, opt.dump_mess);
V_new = dumping(V, V_new, opt.dump_mess);
PRELI = (varNoise + V_new).^(-1);
PRELI2 = (Yeff - W_new) .* PRELI;
var_1 = max(1e-20, sum(PRELI) * Gm2 + HT(PRELI) .* (1 - 2 * Gmean) );
var_2 = HT(PRELI2) - sum(PRELI2) * Gmean;

%% keep track of the previous step
av_mess_old = prior.av_mess;
var_mess_old = prior.var_mess;

%% update
if t > 1 && (strcmp(opt.prior, 'GaussExponential') || strcmp(opt.prior, 'SparseExponential') || strcmp(opt.prior, 'GaussLaplace') )
	computeMeanField;
	prior.R = 1 * (var_2 ./ var_1 + prior.av_mess) + mf * opt.weightMf;
else 
	prior.R = var_2 ./ var_1 + prior.av_mess; 
end
prior.S2 = var_1.^(-1);

%% denoising with prior
prior = F(prior);
V = V_new;
W = W_new;