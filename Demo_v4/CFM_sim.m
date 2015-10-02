function [image, mesure, reconstruction] = CFM_sim(F, object, sparsity, intensity, print)

if nargin < 3, sparsity = 1; end
if nargin < 4, intensity = 2000; end
if nargin < 5, print = 1; end

global n PSF MM filename pathname

switch object
    case 0
        image = double(importdata('real_cam_256.tif'));    
    case 1
        [filename, pathname] = uigetfile({'*.*','All Files'});
        image = double(importdata([pathname, filename]));
    case 2
%% PARAMETERS

n = 512;            % image size
sigma = 1;          % stdev PSF
alpha = 2;
m = n/alpha;

px_size = 800*10^(-3);  % taille du pixel en µm
S = n^2*sparsity/10000;
k = floor(S) ;          % number of beads

%% POINT SPREAD FUNCTION
psf = 2;
switch psf
    case 1
        ksize = 31;
        PSF = zeros(ksize);
        
        s = 1;
        m = ksize/2;
        [X, Y] = meshgrid(1:ksize);
        PSF = (1/(2*pi*s^2))*exp(-((X-m).^2 + (Y-m).^2)/(2*s^2));
    case 2
        PSF = fspecial('gaussian', [256 256], sigma);
end

beads = fspecial('disk', 2); % circular filter

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
        MM = [];                 %  M downsizing (& M'),
        for I = 1:m
            MM = blkdiag(MM,v);
        end;
        
        image = MM*im*MM';
    case 2
        image = imresize(im, 1/alpha, 'bilinear');
end
end

%% Sensing operator (function defined via F = @(z) ... )

if isa(F,'function_handle'),
    b = F(image);
else
    image_vector = reshape(image,[],1);
    b = zeros(size(F,1),1);
    for i = 1:size(F,2)
        proj = image_vector.*F(:,i);
        b(:) = sum(proj);
    end
end

mesure = 1e12*imnoise(b*1e-12, 'poisson');

% [reconstruction] = Fast_IHT_v2(256,mesure,@Afor2f,@Aback2f,1,100,3,1);
% 
% switch print
%     case 1
%         
%     case 2
% %% Plot results
% figure(1); clf;
% 
% subplot(1,3,1);
% imshow(image,[]);  
% if object == 1
%     title('Original image');
% else title ('Simulated image');
% end
% axis square; axis off
% 
% subplot(1,3,2);
% plot(mesure); title('PMT measurements');
% axis square; axis on
% 
% subplot(1,3,3);
% imshow(reconstruction,[]); title('Reconstruction');
% axis square; axis off
% 
% end