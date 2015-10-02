im = double(imread('cells_snap.jpg') );
% im = double(imread('Lena.jpg'));
[M, N] = size(im);

C = fdct_wrapping(im, 1);
percent = 0.1;
Ccut = cell(size(C) );

for j = 1 : numel(C)
	for l = 1 : numel(C{j} )
		[sort_, ind] = sort(abs(C{j}{l}(:) ), 'descend');			
		cut = sort_(ceil(numel(sort_) * percent) );	
		if j > 1
			Ccut{j}{l} = C{j}{l} .* (abs(C{j}{l}) >= cut);
		else
			Ccut{j}{l} = C{j}{l};
		end
	end
end

rec = ifdct_wrapping(Ccut, 1, M, N);

subplot(1, 2, 1); imagesc(im); title('original');
subplot(1, 2, 2); imagesc(rec); title('tresholded');
