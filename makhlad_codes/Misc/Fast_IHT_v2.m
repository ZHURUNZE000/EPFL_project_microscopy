%  function Fast_IHT x = Fast_IHT(N,b,A,At,sp_domain,smooth,positivity,nmax,dlevel,mask)
%
%
%  Input : N - Number of pixels
%          b - measurements (m x 1 vector)
%          A - sensing operator (function defined via A = @(z) ... )
%          At - adjoint sensing operator (function defined via At = @(z) ... )
%
%  Options : nmax - number of iteration
%            Positivity - enforces the positivity of the solution
%            dlevel - detection level - default is 3
%            poisson - set to 1 in case of Poisson noise
%
%  Output : x - N x N reconstructed image
%
%  Example : [x] = Fast_IHT(N,b,A,At,1,500,3,0);
%
%  version : February, 10 2014 - J.Bobin
%


function [x,SigM] = Fast_IHT2(N,b,A,At,positivity,nmax,dlevel,poisson)

%--- Initialization

x = 0.*At(b);
Q = @(z) A(At(z));
L = power_method(Q, numel(b), 25, 1, 0);
alpha = 1/L;

SigM = ones(size(b));

ktf = dlevel;
kt = 15;
dkt = (kt - ktf)/(0.5*nmax);

%--- Main loop

y = x;
xp = x;
tk = 1;

Ay = 0.*b;

for ll=1:nmax 
   
    %--- Compute the gradient
    
    Delta = - alpha*At(SigM.*(b - Ay));
    
    x = y - Delta;
        
    %--- Sparsity constraint
    
    x = x.*(abs(x) > kt*mad(Delta));
                
    %--- Positivity constraint
    
    if positivity x = x.*(x > 0);end
    
    %--- Update
    
    tkp = 0.5*(1+ sqrt(1+4*tk^2));
        
    y = x + (tk-1)/tkp*(x - xp);
    
    xold = y;
    
    xp = x;
    tk = tkp;
    
    kt = max(ktf,kt - dkt);
    
    Ay = A(y);
    
    if poisson
        SigM = abs(Ay) + 7*mad(abs(Ay));
        if max(SigM) == 0
            SigM(:) = 1;
        else
            SigM = 1./SigM;SigM = SigM/max(SigM(:));
        end
    end
            
end

x = reshape(x,N,N);


%############# Auxiliary functions

function [xout,seuil] = Thresh_UDWT(x,N,kt,smooth)

NbrScales = log2(N)-2;
alpha = ones(NbrScales);

seuil = zeros(1,NbrScales);

if smooth 
    alpha = (NbrScales-1:-1:0);
    alpha = alpha/max(alpha);
    alpha = (sqrt(3)-0.1)*alpha+ 0.1;
end

[coeff,xx,N] = UDWT2_ND(reshape(x,N,N),NbrScales, MakeONFilter('Symmlet',4));    
coeff = reshape(coeff,N,[]);
%if nargin > 3 coeff(:,N+1:4*N) = 0;end

g = coeff(:,N+1:end);

for qq = 1:NbrScales 
       
    temp = g(:,3*(qq-1)*N+1:3*qq*N);
    lambda = mad(temp(:));
    seuil(qq) = lambda;
    
    temp(:,1:N) = temp(:,1:N).*(abs(temp(:,1:N)) > kt*alpha(qq)*lambda);
    temp(:,N+1:2*N) = temp(:,N+1:2*N).*(abs(temp(:,N+1:2*N)) > kt*alpha(qq)*lambda);
    temp(:,2*N+1:3*N) = temp(:,2*N+1:3*N).*(abs(temp(:,2*N+1:3*N)) > kt*alpha(qq)*lambda);
    
    g(:,3*(qq-1)*N+1:3*qq*N) = temp;
    
end
coeff(:,N+1:end) = g;
coeff(:,1:N) = coeff(:,1:N).*(abs(coeff(:,1:N)) > 3*mad(reshape(coeff(:,1:N),[],1)));
[xout,img] = IUDWT2_ND(coeff(:),N,N,N,NbrScales,MakeONFilter('Symmlet',4));
xout = xout(:);

%#############

function [coeff,xx,N] = UDWT2_ND(x,scale, qmf);

[n,m] = size(x);
N = n;
M = m;
if (floor(log2(n)) < log2(n)) N = 2^(floor(log2(n))+1); end
if (floor(log2(m)) < log2(m)) M = 2^(floor(log2(m))+1); end

N = max(N,M);

xx = zeros(N,N);
xx(1:n,1:m) = x;

%--- Mirroring

if n > m
   xx(:,m+1:end) = x(:,m:-1:1); 
end

if m > n
   xx(n+1:end,:) = x(n:-1:1,:);
end


[yl,yh,L] = mrdwt(xx,qmf,scale);

coeff = reshape([yl,yh],[],1);


%############# 

function [x,img] = IUDWT2_ND(coeff,N,n,m,scale, qmf);


coeff = reshape(coeff,N,[]);
yl = coeff(:,1:N);
yh = coeff(:,N+1:end);

[img,L] = mirdwt(yl,yh,qmf,scale);

x = img(1:n,1:m);

%#############

% Use power method to estimate maximum eigenvalue of a matrix A
% A could be either a function handle or a real matrix
% n is the size of the matrix
% niter : number of iterations
% nrep is the number of times to run the algorithm.. the result will then
% be the average result over all nrep runs

% Arguments:
%   A - the matrix (either function handle or matrix)
%   n - size of the matrix (A is nxn)
%   niter - How many iterations per run?
%   nrep - How many runs?

function [lambda_max, x ,history] = power_method(A, n, niter, nrep, verbose)

if nargin < 3, niter = 20; end
if nargin < 4, nrep = 1; end
if nargin < 5, verbose = 0; end

history = zeros(1,niter);
estimates = zeros(1,nrep);
lambda_max = 0;
for r=1:nrep
    % pick x randomly
    if verbose,
        fprintf('  power method: repetition %d\n', r);
    end
    x = rand(n,1);
    
    for j=1:niter,
        if isa(A,'function_handle'),
            y = A(x);
        else
            y = A*x;
        end
       x = y / norm(y);
       history(j) = history(j) + norm(y);
    end
    estimates(r) = norm(y);
    lambda_max = lambda_max + norm(y);
end

history = history / nrep;
lambda_max = lambda_max / nrep;
if verbose,
    fprintf('  average estimate of max. eigenvalue = %f\n', mean(estimates));
    fprintf('  power method: standard deviation of estimates = %f\n', std(estimates));
end

