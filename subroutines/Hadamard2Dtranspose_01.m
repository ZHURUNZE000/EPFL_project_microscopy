function y = Hadamard2Dtranspose_01(signal, M, N, IndexRandom)

C = sqrt(N) * (IndexRandom(:, 1) - 1) + IndexRandom(:, 2);
C2 = C(1 : M);
signal2 = zeros(1, N);
signal2(C2) = signal;
image = reshape(signal2 ,sqrt(N), sqrt(N) );
p = hadamards(hadamards(image).');
y = .5 * (sum(signal) + reshape(p, N, 1) ).';

end



