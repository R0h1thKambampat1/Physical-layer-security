M = 4; N = 3; K = 2;
var = 1;
e = ones(M,1);
dB = 0:5:25;
samples= 1e5;
SOP=0;
H = sqrt(var/2)*(randn(N,M)+i*randn(N,M));
G = sqrt(var/2)*(randn(K,M)+i*randn(K,M));
for i = 1:M^2
    