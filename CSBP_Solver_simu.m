function [X, varargout] = CSBP_Solver_simu(Y, H, HT, opt)

initialisation_simu;

%% Construction of the prior dependent class and noise_and_error
prior = Prior(rho, alpha_, opt.learn, opt.prior, opt.dump_learn, R, S2, av_mess, var_mess, param_1, param_2, param_3); 
F = str2func(prior.func);

%% Starting main code
t = 1;
while (t <= opt.tMax)

    AMP_simu;
    
    % Learning of the noise if activated
    if opt.learnNoise
        varNoise = dumping(varNoise, ((Y - W).^2 * (1 + V / varNoise).^(-2).') ./ sum((1 + V / varNoise).^(-1) ), opt.dump_learn);
    end
    
    % print infos and dynamics of reconstruction to screen
    if mod(t, opt.print) == 0        
        conv_ = mean((av_mess_old - prior.av_mess).^2) / mean(prior.av_mess.^2);        
        printToScreen; 
        if (conv_ < opt.conv_)
            pr = sprintf('Converged : convergence = %e', conv_); disp(pr);            
            break;
        end
        if opt.showImage
            imagesc(reshape(prior.av_mess, 256, 256) ); 
            drawnow;
        end
    end

    t = t + 1;
end

X = prior.av_mess;

if strcmp(opt.prior, 'GaussExponential') || strcmp(opt.prior, 'SparseExponential') || strcmp(opt.prior, 'GaussLaplace')
    varargout{1} = prior.weightNoise;
end

end