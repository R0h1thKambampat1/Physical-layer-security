clear all;
close all;
clc;
%% Secure Relay and Jammer Selection for Physical Layer Security

%% Initialisations
% Declaring nodes in the network
M = 4; N = 3; K = 2; % M -> number of intermediate nodes, N -> number of antennas at the destination, K -> is the number of eavesdroppers
SNR_dB = 0:5:25; % No of instances of SNR we would like to plot SNR for
% Generating H matrix - H matrix contains the channel between each intermediate node to ith antenna of the destination.
samples=10^5;
SOP=0;
for ii=1:samples
H = [];
var_mat = []; % a matrix that consists of all the variances used for all the legitimate channels
w = ones(N,1); % initialising beamformer matrix
e = ones(M,1);
gamma_E = [];
gamma_D = [];
var = 1;
G = [];


     G = [];
     H = [];
for j = 1:M
    var = 1;
    var_mat = [var_mat var];
    curr = sqrt(var/2)*(randn(N,1)+i*randn(N,1));
    for k = 1:N
        curr(k,1) = curr(k,1)*curr(k,1)';
    end
    H = [H curr];
end

% Generating G matix - G matrix contains the channel between each intermediate node and ith eavesdropper

for j = 1:M
    curr = sqrt(var/2)*(randn(K,1)+i*randn(K,1));
    for k = 1:K
        curr(k,1) = curr(k,1)*curr(k,1)';
    end
    G = [G curr];
end

gmk = sqrt(1/2)*(randn(M,K)+i*randn(M,K));

for si = 1:length(SNR_dB)
    SNR = 10^(SNR_dB(si)/10);


%     gamma_D_SNR = []; 
%     gamma_E_SNR = []; 
%     for j = 1:2^M
        r_select=randi([1 M],[1 1]);
        b =de2bi(r_select, M);
%         b = de2bi(j-1, M); % b is the matrix that tells us the allocation of jammers and relays among the intermediate nodes
        H_J = H*diag(b); % H_J is the resultant matrix of Jammer nodes
        H_R = H*(eye(M) - diag(b)); % H_R is the resultant matrix of Relay nodes
        temp_gammaD = (SNR*w'*H_R*e*e'*H_R'*w)/(SNR*w'*H_J*H_J'*w +w'*w);
        gamma_D_SNR = temp_gammaD;
        
        
        
        [~,L2] = find(b==1);[~,L1] = find(b==0);
        g_nm = gmk(L1, :);
        g_dm = gmk(L2, :);
        G_NM = abs(sum(g_nm)).^2;
        G_DM = sum(abs(g_dm).^2);
        
        gamma_E_SNR = max(SNR*G_NM./(1+G_DM*SNR));
        
%     end
    gamma_E = [gamma_E gamma_E_SNR];
    gamma_D = [gamma_D gamma_D_SNR];
end

R=max(0.5* (log2(1+gamma_D)-log2(1+gamma_E)),0);
Ri=(R<1);

SOP=SOP+Ri;
end

Pout=sum(SOP)/(samples*2^(M));
    semilogy(SNR_dB, Pout,'-sr','LineWidth',2);
    hold on; grid on
    ylim([10^-4 1])
% R = zeros(2^M-1, SNR_dB);% Secrecy rate
