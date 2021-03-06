function y = Hadamard2D_01(signal, M, N, IndexRandom)

image = reshape(signal, sqrt(N), sqrt(N) );
p = hadamards(hadamards(image).');
y1 = .5 * (sum(signal) + reshape(p, N, 1) );
y = y1(sqrt(N) * (IndexRandom(:, 1) - 1) + IndexRandom(:, 2) );
y = y(1 : M);

end

%Correspondance -> pattern i,j is (i-1)*256+j
%However, things are shifted, so patern 
%index_random_part1(k,1),Index_random_part1(k,2) is
%256*(index_random_part1(k,1)-1)+Index_random_part1(k,2)