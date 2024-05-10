clear all; close all; clc;
%% Generating H and G matrices
M = 4; N = 3; K = 2;
% M is number of intermediate nodes
% N is the number of antennas at the destination
% K is the number of eavesdroppers
% The value of gamma_e cannot be known because of the absence of the CSI of
% eavesdroppers
%%
% samples= 2e4;
% SOP=0;
% for samp = 1:samples
% L2 = 0
var = 1;
jammer_set = [];
relay_set = [1:M];
H = sqrt(var/2)*(randn(N,M)+i*randn(N,M));
G = sqrt(var/2)*(randn(K,M)+i*randn(K,M));
e = ones(M,1);
Q = (H*e)*(H*e)';
[u, v] = eig(Q);
[~, label] = min(abs(eig(Q)));
w = u(:, label);

% L2 = 1
temp = w' * H;
j = min(abs(temp));
[~, label] = find(abs(temp) == j);
jammer_set = [jammer_set, label];

% calculating new w
b = zeros(M, 1);
b(label) = 1;
H_J = H*diag(b);
H_R = H*diag(1-b);
Q = inv(H_J*H_J' + eye(N)) * (H_R*e)*(H_R*e)';
[u, v] = eig(Q);
[~, label] = min(abs(eig(Q)));
w = u(:, label);
jammer_set = [jammer_set label];
relay_set = [1:jammer_set-1, jammer_set+1:M];
% L2 = 2
b = zeros(M, 1);
b(jammer_set) = 1;
H_J = H*diag(b);
H_R = H*diag(1-b);
Q = inv(H_J*H_J' + eye(N)) * (H_R*e)*(H_R*e)';
[u, v] = eig(Q);
[~, label] = min(abs(eig(Q)));
w = u(:, label);



    
    
    