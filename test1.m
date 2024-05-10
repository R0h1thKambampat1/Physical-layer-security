clear all;
close all;
clc;
%% Secure Relay and Jammer Selection for Physical Layer Security

%% Initialisations
% Declaring nodes in the network
M = 4; N = 3; K = 2; % M -> number of intermediate nodes, N -> number of antennas at the destination, K -> is the number of eavesdroppers
SNR_upperlimit = 25; % No of instances of SNR we would like to plot SNR for
% Generating H matrix - H matrix contains the channel between each intermediate node to ith antenna of the destination.
H = [];
var_mat = []; % a matrix that consists of all the variances used for all the legitimate channels
for j = 1:M
    var = 0.5;
    var_mat = [var_mat var];
    curr = sqrt(var/2)*(randn(N,1)+i*randn(N,1));
    for k = 1:N
        curr(k,1) = curr(k,1)*curr(k,1)';
    end
    H = [H curr];
end

% Generating G matix - G matrix contains the channel between each intermediate node and ith eavesdropper
var = 0.1;
G = [];
for j = 1:M
    curr = sqrt(var/2)*(randn(K,1)+i*randn(K,1));
    for k = 1:K
        curr(k,1) = curr(k,1)*curr(k,1)';
    end
    G = [G curr];
end

%% Calculating SNRs at Destination 
w = ones(N,1); % initialising beamformer matrix
e = ones(M,1);

gamma_D = [];
for SNR = 0:SNR_upperlimit-1
    gamma_D_SNR = []; 
    for j = 0:2^M-1
        b = de2bi(j, M); % b is the matrix that tells us the allocation of jammers and relays among the intermediate nodes
        H_J = H*diag(b); % H_J is the resultant matrix of Jammer nodes
        H_R = H*(eye(M) - diag(b)); % H_R is the resultant matrix of Relay nodes
        temp_gammaD = (SNR*w'*H_R*e*e'*H_R'*w)/(SNR*w'*H_J*H_J'*w +w'*w);
        gamma_D_SNR = [gamma_D_SNR; temp_gammaD];
    end
    gamma_D = [gamma_D gamma_D_SNR];
end
% gamma_D is 2^M x SNR_upper limit with (i,j) element representing gamma_D
% for ith configuration and jth SNR value that is set.

%% Rate at destination
Rd = [];
for i = 1:2^M
    for j = 1:SNR_upperlimit
        Rd(i,j) = 0.5*log2(1+gamma_D(i,j));
    end
end
for i = 1:2^M
    figure;
    plot([0:24], Rd(i,:));
    title (sprintf('%dth configuration',i))
    
end
%% Calculating SNRs at Eavesdroppers

% gamma_E = [];
% for SNR = 0:SNR_upperlimit-1
% gamma_E_temp = [];
% for config = 0:2^M-1
%     b = de2bi(config, M);
%     gam = 0;
%     for k = 1:K
%         denom = 0; num = 0;
%         for j = 1:M
%             if(b(1, M) == 1)
%                 denom = denom+G(k, j)*G(k,j)';
%             else 
%                 num = num+G(k,j);
%             end
%         end
%         temp = SNR*(num*num')/(SNR*denom + 1);
%         gam = max(gam, temp);% the gamma value for each configuration, temp is the value of gamma at kth eavesdropper. A max of gammas at all eavesdropper for a certain configuration is considered
%     end
% gamma_E_temp = [gamma_E_temp; gam]; % gamma_E_temp is the column vector of gam s for each configuration.
% end
% gamma_E = [gamma_E gamma_E_temp];
% end
% 
% % SNR = 25 ;
% % 
% % gamma_temp =[];
% % for j = 0:2^M-1
% % b = de2bi(j, M);
% % gamma_max = 0;
% % for k = 1:K
% %     num=0; den =0;
% % for i = 1:M
% %     
% %     if(b(1, M) == 1)
% %         den = den + G(i,k)*G(i,k)';
% %     else 
% %         num = num + G(i,k);
% %     end
% % end
% % gamma_eavs = (SNR*num*num')/(SNR*den+1);
% % gamma_max = max(gamma_max, gamma_eavs);
% % end
% % gamma_temp = [gamma_temp; gamma_max]
% % end
% 
% 
% %% Calculating Secrecy Rate
% R = zeros(2^M, SNR_upperlimit);% Secrecy rate
% for j = 1:2^M
%     for k = 1:SNR_upperlimit
%         R(j,k) = max(0.5*log2((1+gamma_D(j,k))/(1+gamma_E(j,k))), 0);
%     end
% end
% %% SOP calculation
% sop = [];
% 
% for k = 1:SNR_upperlimit
%     count = 0;
%     for j = 1:2^M
%         if(R(j,k) < 10^0.3 -1)
%             count = count+1;
%         end
%     end
%     sop = [sop count/(2^M)];
% end
% 







