classdef Prior
    % This class contains all the prior-dependent functions including learnings
    
    properties
        av_mess;
        var_mess; 
        R; 
        S2; 
        rho; 
        learn;         
        alpha; 
        dump_learn; 
        param_1;
        param_2;
        param_3;
        func;
        weightNoise;
    end
    
    methods
        
        function prior = Prior(rho, alpha, learn, choice_prior, dump_learn, R, S2, av_mess, var_mess, param_1, param_2, param_3)
            % Constructor function
            prior.R = R; 
            prior.S2 = S2; 
            prior.rho = rho; 
            prior.learn = learn; 
            prior.alpha = alpha; 
            prior.av_mess = av_mess;             
            prior.var_mess = var_mess;             
            prior.dump_learn = dump_learn; 

            switch choice_prior                             
                case 'L1'
                    prior.param_1 = min(param_1, param_2);
                    prior.param_2 = max(param_1, param_2);
                    prior.func = 'PriorL1';
                    disp('L1')      
                case 'GaussExponential'
                    prior.param_1 = param_1;
                    prior.param_2 = param_2;
                    prior.param_3 = param_3;
                    prior.func = 'PriorGE';
                    disp('GaussExponential')
                case 'SparseExponential'
                    prior.param_1 = param_1;
                    prior.param_2 = param_2;
                    prior.param_3 = param_3;
                    prior.func = 'PriorSE';
                    disp('SparseExponential')
                case 'GaussLaplace'
                    prior.param_1 = param_1;
                    prior.param_2 = param_2;
                    prior.param_3 = param_3;
                    prior.func = 'PriorGLap';
                    disp('GaussLaplace')                          
            end
        end                              

        function prior = PriorGE(prior)
            % Gauss-Exponential prior : p(x) ~ exp(-(x - mG)^2 / (2 * varG) ) / sqrt(2 * pi * varG) + rho * I(x > 0) * expo * exp(-expo * x), expo > 0 : param_1 = expo; param_2 = mG; param_3 = varG;
            R_ = prior.R; 
            S2_ = prior.S2;
            rho_ = prior.rho; 
            expo_ = prior.param_1; 
            mG = prior.param_2; 
            varG = prior.param_3;

            a = exp(-R_.^2 ./ (2 .* S2_) );
            expo1 = exp(-.5 .* (R_ - mG).^2 ./ (varG + S2_) );
            f1 = expo1 ./ sqrt(2 .* pi .* (varG + S2_) );                        
            b = exp(-expo_ .* R_ + expo_.^2 .* S2_ ./ 2);
            c = erfc((expo_ .* sqrt(S2_) - R_ ./ sqrt(S2_) ) ./ sqrt(2) );
            Z = f1 + rho_ .* expo_ .* b ./ 2 .* c;            
            h1 = expo1 .* (mG.^2 .* S2_.^2 + S2_ .* varG .* (2 .* mG .* R_ + S2_) + (R_.^2 + S2_) .* varG.^2) ./ (sqrt(2 .* pi) .* (varG + S2_).^(5 ./ 2) );
            h2 = S2_ ./ sqrt(2 .* pi) .* (-expo_ .* sqrt(S2_) + R_ ./ sqrt(S2_) ) .* a + b ./ 2 .* c .* (S2_ + (R_ - expo_ .* S2_).^2);            
            prior.av_mess = ((mG .* S2_ + R_ .* varG) ./ (varG + S2_) .* f1 + rho_ .* expo_ .* (sqrt(S2_ ./ (2 .* pi) ) .* a + (R_ - expo_ .* S2_) ./ 2 .* b .* c) ) ./ Z;
            prior.var_mess = max(1e-20, (h1 + rho_ .* expo_ .* h2) ./ Z - prior.av_mess.^2);                        
            prior.weightNoise = f1 ./ Z;

            if prior.learn     
                % Expectation maximisation                                           
                prior.rho = dumping(prior.rho, mean(prior.weightNoise < .5), prior.dump_learn);
                % Expectation maximisation
                prior.param_1 = dumping(prior.param_1, max(1e-10, 1 / mean(prior.av_mess(prior.weightNoise < .5) ) - 0 ), prior.dump_learn);
                % Expectation maximisation
                prior.param_2 = dumping(prior.param_2, mean(prior.av_mess(prior.weightNoise > .5) ), prior.dump_learn);
                % % Expectation maximisation
                % prior.param_3 = dumping(prior.param_3, mean(prior.var_mess(prior.weightNoise > .5) ), prior.dump_learn);

                % % Bayes optimal
                % prior.param_2 = dumping(prior.param_2, sum(prior.weightNoise .* R_ ./ (S2_ + varG) ) / sum(prior.weightNoise ./ (S2_ + varG) ), prior.dump_learn);
                % % Bayes optimal                
                % prior.param_3 = dumping(prior.param_3, sum(prior.weightNoise .* ((mG - R_) ./ (varG + S2_) ).^2) / sum(prior.weightNoise .* varG ./ (varG + S2_) ), prior.dump_learn);                
            end                      
        end 

        function prior = PriorSE(prior)
            % Exponential sparse prior : p(x) ~ (1 - rho) * delta(x) + rho * I(x > 0) * expo * exp(-expo * x), expo > 0 : param_1 = expo;
            R_ = prior.R; 
            S2_ = prior.S2; 
            rho_ = prior.rho; 
            expo_ = prior.param_1;
            
            a = exp(-R_.^2 ./ (2 .* S2_) );
            b = exp(-expo_ .* R_ + expo_.^2 .* S2_ ./ 2);
            c = erfc((expo_ .* sqrt(S2_) - R_ ./ sqrt(S2_) ) ./ sqrt(2) );
            d = (1 - rho_) ./ sqrt(2 .* pi .* S2_) .* a;
            Z = d + rho_ .* expo_ .* b ./ 2 .* c;
            prior.av_mess = rho_ .* expo_ .* (sqrt(S2_ ./ (2 .* pi) ) .* a + (R_ - expo_ .* S2_) ./ 2 .* b .* c) ./ Z;            
            prior.var_mess = max(1e-20, rho_ .* expo_ .* (S2_ ./ sqrt(2 .* pi) .* (-expo_ .* sqrt(S2_) + R_ ./ sqrt(S2_) ) .* a + b ./ 2 .* c .* (S2_ + (R_ - expo_ .* S2_).^2 ) ) ./ Z - prior.av_mess.^2);
            prior.weightNoise = d ./ Z;

            if prior.learn                
                % Expectation maximisation
                prior.param_1 = dumping(prior.param_1, max(1e-10, 1 / mean(prior.av_mess(prior.weightNoise < .5) ) ), prior.dump_learn);                
                % Expectation maximisation                                           
                prior.rho = dumping(prior.rho, mean(prior.weightNoise < .5), prior.dump_learn); 

                % % Bayes optimal
                % Z_rho = (1 - rho_) .* a ./ sqrt(2 .* pi .* S2_) + (prior.av_mess - rho_ .* sqrt(S2_ ./ (2 .* pi) .* a) ) ./ (R_ - expo_ .* S2_);
                % prior.rho = dumping(prior.rho, abs(((prior.av_mess - rho_ .* sqrt(S2_ ./ (2 .* pi) .* a) ) ./ (R_ - expo_ .* S2_) * Z_rho.^(-1).') ./ (a ./ sqrt(2 .* pi .* S2_) * Z_rho.^(-1).') ), prior.dump_learn);               
            end            
        end

        function prior = PriorGLap(prior)
            % Gauss-Laplace prior : p(x) ~ exp(-(x - mG)^2 / (2 * varG) ) / sqrt(2 * pi * varG) + rho * expo / 2 * exp(-expo * |x|), expo > 0 : param_1 = expo; param_2 = mG; param_3 = varG;
            R_ = prior.R; 
            S2_ = prior.S2;
            rho_ = prior.rho; 
            expo_ = prior.param_1; 
            mG = prior.param_2; 
            varG = prior.param_3;

            n = exp(-.5 .* (mG - R_).^2 ./ (varG + S2_) );
            exp1m = exp(.5 .* expo_ .* (expo_ .* S2_ - 2 .* R_) );
            exp1p = exp(.5 .* expo_ .* (expo_ .* S2_ + 2 .* R_) );
            erfcm = erfc((expo_ .* S2_ - R_) ./ sqrt(2 .* S2_) );
            erfcp = erfc((expo_ .* S2_ + R_) ./ sqrt(2 .* S2_) );
            exp0 = exp(-.5 .* R_.^2 ./ S2_);
            Z = n ./ sqrt(2 .* pi .* (varG + S2_) ) + rho_ ./ 4 .* expo_ .* (exp1m .* erfcm + exp1p .* erfcp); 
            prior.av_mess = ((R_ .* varG + mG .* S2_) ./ (varG + S2_) .* n ./ sqrt(2 .* pi .* (varG + S2_) ) + rho_ ./ 4  .* expo_ .* (exp1m .* erfcm .* (R_ - expo_ .* S2_) + exp1p .* erfcp .* (R_ + expo_ .* S2_) ) ) ./ Z;
            fb = (n ./ sqrt(2 .* pi .* (varG + S2_).^5) .* ((mG .* S2_).^2 + S2_ .* varG .* (2 .* mG .* R_ + S2_) + varG.^2 .* (R_.^2 + S2_) ) + rho_ ./ 4 .* expo_ .* (-2 .* exp0 .* expo_ .* S2_ .* sqrt(2 .* S2_ ./ pi) + exp1m .* erfcm .* (S2_ + (R_ - S2_ .* expo_).^2) + exp1p .* erfcp .* (S2_ + (R_ + S2_ .* expo_).^2) ) ) ./ Z;                   
            prior.var_mess = max(1e-20, fb - prior.av_mess.^2);                        
            prior.weightNoise = n ./ sqrt(2 .* pi .* (varG + S2_) ) ./ Z;

            if prior.learn     
                % Expectation maximisation                                           
                prior.rho = dumping(prior.rho, mean(prior.weightNoise < .5), prior.dump_learn);
                % Expectation maximisation
                % prior.param_1 = dumping(prior.param_1, max(1e-10, 1 / mean(prior.av_mess(prior.weightNoise < .5) ) ), prior.dump_learn);
                % Expectation maximisation
                prior.param_2 = dumping(prior.param_2, mean(prior.av_mess(prior.weightNoise > .5) ), prior.dump_learn);
                % % Expectation maximisation
                % prior.param_3 = dumping(prior.param_3, mean(prior.var_mess(prior.weightNoise > .5) ), prior.dump_learn);
            end                      
        end

        function prior = PriorL1(prior)
            % L1 optimization (soft tresholding) : p(x) ~ lim_{beta -> infinity} exp{-beta * |x|}, where the x values are bounded by [min, max] : param_1 = min; param_2 = max;
            min_ = prior.param_1; 
            max_ = prior.param_2;
            R_ = prior.R; 
            S2_ = prior.S2;

            prior.av_mess = min(max_, (R_ > 0) .* (R_ - S2_) .* (R_ > S2_) ) + max(min_, (R_ < 0) .* (R_ + S2_) .* (-R_ > S2_) );
            prior.var_mess = S2_ .* (abs(R_) > S2_);            
        end                
        
    end
    
end