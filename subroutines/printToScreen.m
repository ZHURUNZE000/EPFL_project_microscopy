switch opt.prior    
    case ('SparseExponential')
        disp('t   noise   expo   rho   conv');         
        disp(sprintf('%d %.2e %.2e %.2e %.2e %.2e %.2e', [t, varNoise, prior.param_1, prior.rho, conv_] ) );        
    case ('GaussExponential')
        disp('t   noise    expo      mG      vG      rho      conv');         
        disp(sprintf('%d %.2e %.2e %.2e %.2e %.2e %.2e', [t, varNoise, prior.param_1, prior.param_2, prior.param_3, prior.rho, conv_] ) );
    case ('GaussLaplace')
        disp('t   noise   expo   mG   vG   rho   conv');         
        disp(sprintf('%d %.2e %.2e %.2e %.2e %.2e %.2e', [t, varNoise, prior.param_1, prior.param_2, prior.param_3, prior.rho, conv_] ) );        
    case ('L1')
        disp('t   noise   rho   conv');         
        disp(sprintf('%d %.2e %.2e %.2e %.2e %.2e %.2e', [t, varNoise, prior.rho, conv_,] ) );        
end