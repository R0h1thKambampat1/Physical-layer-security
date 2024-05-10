clear all; close all; clc;
%%
M = 4; N = 3; K = 2;
var = 1;
e = ones(M,1);
dB = 0:5:25;
samples= 1e4;
SOP=0;
b = [0 0 0 1; 0 0 1 0; 0 0 1 1; 0 1 0 0;
    0 1 0 1; 0 1 1 0; 0 1 1 1; 1 0 0 0;
    1 0 0 1; 1 0 1 0; 1 0 1 1; 1 1 0 0;
    1 1 0 1; 1 1 1 0; 1 1 1 1];
R_all_samps = [];
for iii = 1:samples
    H = sqrt(var/2)*(randn(N,M)+i*randn(N,M));
    G = sqrt(var/2)*(randn(K,M)+i*randn(K,M));
    
    R_all_SNRs = [];
    for ii = 1:length(dB)
        SNR = 10^(dB(ii)/10);
        R = 0;
        for i = 1:length(b)-1
            jammers = diag(b(i, :));
            relays = eye(M)-jammers;
            H_R = H*relays;
            H_J = H*jammers;
            Q = inv(SNR*(H_J*H_J') + eye(N)) * SNR*(H_R*e)*(H_R*e)';
            gd_temp = max(abs(eig(Q)));
            [~,L2] = find(jammers==1);[~,L1] = find(relays==1);
            ge_num = G(:, L1);
            ge_denom = G(:, L2);
            ge_num = abs(sum(ge_num)).^2;
            ge_denom = sum(abs(ge_denom).^2);
            ge_temp = max(SNR*ge_num./(1+ge_denom*SNR));
            temp = max(0.5* (log2(1+gd_temp)-log2(1+ge_temp)),0);
            R = max(R, temp);
        end
        R_all_SNRs = [R_all_SNRs; R];
    end
    R_all_samps = [R_all_samps R_all_SNRs];
end
R_res = (R_all_samps<1);
R_res = sum(R_res');
R_res = R_res/samples;
semilogy(dB, R_res);