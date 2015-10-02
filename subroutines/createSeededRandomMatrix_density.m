function [G] = createSeededRandomMatrix_density(J, density, Mblock, Nblock)
% Creates a random seeded real binary matrix from the J matrix.

[numBlockL, numBlockC] = size(J);
N = Nblock .* numBlockC;

G = randn(sum(Mblock), N) > 0;
GJ = full(sprandn(sum(Mblock), N, density) ) ~= 0;

% "seeding" of the matrix
startL = 1; stopL = Mblock(1);
for l = 1 : numBlockL
    for c = 1 : numBlockC; 

    	if (J(l, c) == 0)
    		G(startL : stopL, Nblock .* (c - 1) + 1 : Nblock .* c) = 0;
    	elseif (J(l, c) < 1) && (J(l, c) > 0)
			G(startL : stopL, Nblock .* (c - 1) + 1 : Nblock .* c) = GJ(startL : stopL, Nblock .* (c - 1) + 1 : Nblock .* c);
		end    	    

    end

	if (l < numBlockL); startL = stopL + 1; stopL = stopL + Mblock(l + 1); end
    if (l == numBlockL); startL = stopL + 1; stopL = sum(Mblock); end

end

end