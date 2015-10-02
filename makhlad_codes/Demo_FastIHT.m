addpath('Misc');
load acquisition.mat
b = data2b';

[x_FastIHT_8192] = Fast_IHT_v2(256,b,@Afor2f_8192,@Aback2f_8192,1,500,3,1);

[x_FastIHT_4096] = Fast_IHT_v2(256,b(1:4096),@Afor2f_4096,@Aback2f_4096,1,500,3,1);

[x_FastIHT_2048] = Fast_IHT_v2(256,b(1:2048),@Afor2f_2048,@Aback2f_2048,1,500,3,1);

[x_FastIHT_1024] = Fast_IHT_v2(256,b(1:1024),@Afor2f_1024,@Aback2f_1024,1,500,3,1);

[x_FastIHT_512]  = Fast_IHT_v2(256,b(1:512), @Afor2f_512,@Aback2f_512,1,500,3,1);