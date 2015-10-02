%% SENSING MATRIX F

N = 256^2;
M = floor(.1*N);

F = rand(N, M);    

[image, mesure] = CFM_sim(F, 0);