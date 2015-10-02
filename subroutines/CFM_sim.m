function [image, mesure] = CFM_sim(F, sparsity, intensity)

global n PSF M

%% PARAMETERS

n = 512;            % image size
sigma = 1;          % stdev PSF
alpha = 2;
m = n/alpha;

px_size = 800*10^(-3);  % taille du pixel en µm
S = n^2*sparsity;
k = floor(S) ;          % number of beads

%% POINT SPREAD FUNCTION
psf = 1;

switch psf
    case 1
        %Blur Kernel
        ksize = 31;
        PSF = zeros(ksize);
        
        %Gaussian Blur
        s = 1;
        m = ksize/2;
        [X, Y] = meshgrid(1:ksize);
        PSF = (1/(2*pi*s^2))*exp(-((X-m).^2 + (Y-m).^2)/(2*s^2));
        
    case 2
        PSF = fspecial('gaussian', [256 256], sigma);
end

beads = fspecial('disk', 2);

%% VIRTUAL OBJECT IMAGE

B = zeros(n,n);
support = randsample(n*n,k);
B(support) = randsample(intensity-(intensity*0.1):0.01:intensity+(intensity*0.1),k);
im = imfilter(B,beads);
im = imfilter(im,PSF);

resize = 2;

switch resize
    case 1
        %% RESIZE MATRICE n*n to m*m
        v = ones(1,alpha);
        M = [];                 %  M downsizing (& M'),
        for I = 1:m
            M = blkdiag(M,v);
        end;
        
        image = M*im*M';
    case 2
        image = imresize(im, 1/alpha, 'bilinear');
end

image = reshape(image, [], 1);

%% Sensing operator (function defined via F = @(z) ... )

if isa(F,'function_handle'),
    b = F(image);
else
    b = F*image;
end

mesure = 1e12*imnoise(b*1e-12, 'poisson');
