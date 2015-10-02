function [image, mesure] = CFM_sim_v2(F, object, sparsity, intensity, print)

if nargin < 3, sparsity = 1; end
if nargin < 4, intensity = 2000; end
if nargin < 5, print = 1; end

global n PSF M filename pathname

switch object
    case 1
        % [filename, pathname] = uigetfile({'*.*','All Files'});
        % image = double(importdata([pathname, filename]));
        image = double(importdata('~/Documents/Research/CS/bordeaux_new_mars2015/demos/Demov2/real_cam_256.tif'));
    case 2
%% PARAMETERS

    n = 512;            % image size
    sigma = 1;          % stdev PSF
    alpha = 2;
    m = n/alpha;

    px_size = 800*10^(-3);  % taille du pixel en µm
    S = n^2*sparsity;
    k = floor(S)        % number of beads

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
            M = [];                 %  M downsizing (& M'),
            for I = 1:m
                M = blkdiag(M,v);
            end;
            
            image = M*im*M';
        case 2
            image = imresize(im, 1/alpha, 'bilinear');
    end
end

%% Sensing operator (function defined via F = @(z) ... )

if isa(F,'function_handle'),
    b = F(image(:));
else
%     image = reshape(b,1,[]);
%     b = image*F;
end

plot(b)
drawnow
pause
% mesure = 1e12*imnoise(b*1e-12, 'poisson');
mesure = imnoise(uint8(b), 'poisson');
plot(mesure)
drawnow
pause
mesure = double(mesure);

% [reconstruction] = Fast_IHT_v2(256,mesure,@Afor2f,@Aback2f,1,100,3,1);

% switch print
%     case 1
        
%     case 2
% %% Plot results
% figure(1); clf;

% subplot(1,3,1);
% imshow(image,[]);  
% if object == 1
%     title('Original image');
% else title ('Simulated image');
% end
% axis square; axis off

% subplot(1,3,2);
% plot(mesure); title('PMT measurements');
% axis square; axis on

% subplot(1,3,3);
% imshow(reconstruction,[]); title('Reconstruction');
% axis square; axis off

% end