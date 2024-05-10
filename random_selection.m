clear all; close all; clc;
%% Generating H and G matrices
M = 4; N = 3; K = 2;
% M is number of intermediate nodes
% N is the number of antennas at the destination
% K is the number of eavesdroppers
% The value of gamma_e cannot be known because of the absence of the CSI of
% eavesdroppers
samples= 1e3;
SOP=0;
for samp = 1:samples
    var = 1;
    H = sqrt(var/2)*(randn(N,M)+i*randn(N,M));
    G = sqrt(var/2)*(randn(K,M)+i*randn(K,M));
    %% Finding gamma_D_max and w
    dB = 0:5:25;
    gamma_D_max = [];
    w_res = [];
    for ii = 1:length(dB)
        SNR = 10^(db(ii)/10);
        temp = [];
        for i = 1:2^M-1
            b = de2bi(i-1, M);
            H_R = H*diag(1-b);
            H_J = H*diag(b);
            Q = inv(SNR*H*diag(b)*H' + eye(N)) * (SNR*H*diag(1-b)*H');
            [u, v] = eig(Q);
            [l,m] = find(v==max(max(eig(Q))));
            w{i} = u(:,l);
            temp = [temp abs(max(eig(Q)))];
        end
        gamma_D_max = [gamma_D_max; temp];
        w_res = [w_res; w];
    end
    % Sometimes getting w_res{i,1} as []. To deal with this I am modifying all
    % the elements in the first column (all relays case for all dB) as the same
    % result
    for t = 1:length(dB)
        if(~isempty(w_res{t, 1}))
            break;
        end
    end
    for t1 = 1:length(dB)
        w_res{t1,1} = w_res{t, 1};
    end
    %% Finding gamma_D and gamma_E
    % gamma_D = [];
%     gamma_D = gamma_D_max(:,10);
    gamma_E = [];e = ones(M,1);
    for si = 1:length(dB)
        SNR = 10^(dB(si)/10);
        r_select = randi([1 14], [1 1]);
        b = de2bi(r_select, M);
        %     H_J = H*diag(b);
        %     H_R = H*(eye(M) - diag(b)); % H_R is the resultant matrix of Relay nodes
        %     temp_gammaD = (SNR*w'*H_R*e*e'*H_R'*w_res)/(SNR*w'*H_J*H_J'*w_res +w'*w);
        %     gamma_D_SNR = temp_gammaD;
        
        
        [~,L2] = find(b==1);[~,L1] = find(b==0);
        G_num = G(:, L1);
        G_denom = G(:, L2);
        G_num = abs(sum(G_num)).^2;
        G_denom = sum(abs(G_denom).^2);
        
        gamma_E_SNR = max(SNR*G_num./(1+G_denom*SNR));
        gamma_E = [gamma_E; abs(gamma_E_SNR)];
    end
    gamma_D = gamma_D_max(:, r_select);
    
    R=max(0.5* (log2(1+gamma_D)-log2(1+gamma_E)),0);
    Ri=(R<0.5);
    SOP = SOP+Ri;
end
%%
Pout=SOP/(samples);
semilogy(dB, Pout,'-sr','LineWidth',4);
hold on; grid on;
ylim([10^-4 1]);
xlabel("SNR in dB");
ylabel("SOP");
