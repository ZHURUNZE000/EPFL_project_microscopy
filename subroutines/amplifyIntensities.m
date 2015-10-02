function X = amplifyIntensities(Y, n)

X = Y;

vec = find(X ~= 0);
val = .5;
X(vec) = 1;
for i = 1 : n
	X(min(vec + i, 256^2) ) = val / i;
	X(max(vec - i, 1) ) = val / i;
	X(max(vec - 256 * i, 1) ) = val / i;
	X(min(vec + 256 * i, 256^2) ) = val / i;
end
X(vec) = 1;

end