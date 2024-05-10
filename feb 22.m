clear all; clc;
N = 3; % number of antennas at the destination
M = 4; % number of intermediate nodes
K = 2; % number of eavesdroppers 
var = 2;
P = 10;
H = [];

for i = 1:M
    temp = sqrt(var/2)*(randn(N,1)+i*randn(N,1));
    H = [H temp];
end
G = [];
for i = 1:M
    temp = sqrt(var/2)*(randn(K,1)+i*randn(K,1));
    G = [G temp];
end
w = ones(N,1);
e = ones(M,1);
P = 10; No = 2;

for i = 1:2^M-1
    count = 0;
    temp = de2bi(i, M);
    H_J = H*diag(temp);
    H_R = H*(eye(M) - diag(temp));
    G_J = G*diag(temp);
    G_R = H*(eye(M) - diag(temp));
    gamma_D = (P*w'*H_R*e*e'*H_R'*w)/(P*w'*H_J*H_J'*w + No*w'*w);
    
end
gamma_E = (P*G_R*G_R')/(P*G_J*G_J' + No);
    rate = 0.5*log((1+gamma_D)/(1+gamma_E));
    if (rate<=0)
        count = count + 1;
    end
    
for i = 1:2^M-1
    temp = de2bi(i, M);
   y

