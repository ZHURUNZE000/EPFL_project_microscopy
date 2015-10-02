function [psnr, solve_time, tval3_image] = tval3_experiment(y,A,A_T,original_image,mu,firstMode)
    % Recalculating because I don't like passing a million parameters
    Xmin = min(original_image(:));
    Xmax = max(original_image(:));
    % [num_rows num_cols] = size(original_image);
    num_rows = 256;
    num_cols = 256;
    
    % Put together the function handle structure
    B = @(z,t) PROJ(z,t,A,A_T);
    
    % Set options    
    opts.mu = mu;
    opts.beta = 2^10;
    opts.tol = 1E-6;
    opts.maxit = 1000;
    opts.TVnorm = 1;
    if firstMode; opts.nonneg = true; else opts.nonneg = false; end 
    opts.init = reshape(A_T(y), 256, 256);
    opts.isreal = true;
    

    % Calculate recovery
    tic
    [tval3_image, out] = TVAL3(B,y,num_rows,num_cols,opts);
    solve_time = toc;
    % psnr = PSNR(original_image(:),(Xmax-Xmin)*tval3_image(:) + Xmin);
    psnr = 1;
    
    
function out  = PROJ(in,type,func,func_T)
    if(type == 1)        
        out = func(in);
    else
        if(type == 2)
            out = func_T(in);
        end
    end