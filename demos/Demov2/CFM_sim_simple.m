function [im, meas, noiselessMeas, noise] = CFM_sim_simple(F, sparsity, intensity, sizeBeads, noise_)

n = 512;
k = floor(n^2 * sparsity); % number of beads

PSF = fspecial('gaussian', [256, 256], 1); % point spread function

im = zeros(n);
supp = randsample(n^2, k); % select the random positions of the beads

im(supp) = randsample([intensity * .9 : 0.001 : intensity * 1.1], k); % put the beads on the plane

beads = fspecial('disk', sizeBeads); % filter to make the beads as disks
im = imfilter(im, beads);

im = imfilter(im, PSF); % apply the Gaussian PSF with std = 1

im = imresize(im, .5, 'bilinear'); % shrinks the image by a factor 2

noiselessMeas = F(im(:).'); % noiseless measurement

if noise_
    meas = 1e12 * imnoise(noiselessMeas * 1e-12, 'poisson'); % poisson noise
else
    meas = noiselessMeas;
end

noise = meas - noiselessMeas;

end