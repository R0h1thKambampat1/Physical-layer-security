clear all; close all; clc;
%%
M = 4; N = 3; K = 2; % M is number of intermediate nodes, N is the number of antennas at the destination, K is the number of eavesdroppers
samples = 10^2;
dB = 0:5:25;
SOP=0;
for sample = 1:samples
    var = 1;
    H = sqrt(var/2)*(randn(N,M)+i*randn(N,M));
    G = sqrt(var/2)*(randn(K,M)+i*randn(K,M));
    
    w = [];
    gd = [];
    ge = [];
    
    for i = 1:length(dB)
        % Finding Gamma_D_max(gd), w vectors and ge
        SNR = 10^(dB(i)/10);
        gd_temp = [];
        ge_temp = [];
        for j = 2:15
            % gd calculation
            b = de2bi(j-1, M);
            H_R = H*diag(1-b);
            H_J = H*diag(b);
            Q = inv(SNR * (H_R * H_R') + eye(N)) *  (H_J * H_J');
            [u,v] = eig(Q);
            v_max = abs(max(max(v)));
            [a,~] = (find(abs(v) == v_max));
            curr{j-1} = u(:, a); % current eigen vector corresponding to the max eigen value
            gd_temp = [gd_temp v_max];
            % ge calculation
            [~, L2] = find(b==1); % identifying jammer nodes
            [~, L1] = find(b==0); % identifying relay nodes
            G_num = G(:, L1);
            G_denom = G(:, L2);
            G_num = abs(sum(G_num)).^2;
            G_denom = sum(abs(G_denom).^2);
            temp = max(SNR*G_num./(1+G_denom*SNR));
            ge_temp = [ge_temp abs(temp)];
        end
        gd = [gd; gd_temp];
        ge = [ge; ge_temp];
        w = [w; curr];
    end
    
    for i = 1:14
        R=max(0.5*(log2(1+gd(:,i))-log2(1+ge(:,i))),0);
        Ri=(R<1);
        SOP = SOP+Ri;
    end
end
%%
Pout = SOP/(14*samples);
semilogy(dB, Pout,'-sr','LineWidth',4);
hold on; grid on;
ylim([10^-4 1]);
xlabel("SNR in dB");
ylabel("SOP");

