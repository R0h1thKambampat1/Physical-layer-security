clear all;
close all;
clc;
%%
M = 4; N = 3; K = 2;
dB = 0:5:25;
samples = 10^5;
SOP = [];
% for ii = 1:samples
e = ones(M,1);
hnm = sqrt(1/2)*(randn(N,M)+i*randn(N,M));
gamma_D = [];
gmk = sqrt(1/2)*(randn(M,K)+i*randn(M,K));
for si = 1:length(dB)
    SNR = 10^(dB(si)/10);
    temp = [];
    for i = 1:2^M-1
        b = de2bi(i-1, M);
        Q = inv(SNR*hnm*diag(b)*hnm' + eye(N)) * (SNR*hnm*((1-b)*(1-b)')*hnm');
        temp = [temp; max(eig(Q))];
    end
    gamma_D = [gamma_D, temp];
    
end
        
    
