clear all; close all; clc;
%%
M = 4; N = 3; K = 2;
var = 1;
e = ones(M,1);
dB = 0:5:25;
samples= 1e3;
SOP=0;
for ii =  1:samples
    gd = [];
    ge = [];
    for si = 1:length(dB)
        SNR = 10^(dB(si)/10);
        gd_temp = [];
        ge_temp = [];
        var = 1;
        H = sqrt(var/2)*(randn(N,M)+i*randn(N,M));
        G = sqrt(var/2)*(randn(K,M)+i*randn(K,M));
        w = [];
        for i = 1:2^M -2
            jammer_set = de2bi(i, M);
            relay_set = 1 - jammer_set;
            H_R = H * diag(relay_set);
            H_J = H * diag(jammer_set);
            Q = inv(SNR*(H_J * H_J') + eye(N)) * (SNR*(H_R*e)*(H_R*e)');
            [u, v] = eig(Q);
            [u,v] = eig(Q);
            v_max = abs(max(max(v)));
            [a,~] = (find(abs(v) == v_max));
            gd_temp = [gd_temp v_max];
            
            [~,L2] = find(jammer_set==1);
            [~,L1] = find(relay_set==1);
            ge_num = G(:, L1);
            ge_denom = G(:, L2);
            ge_num = abs(sum(ge_num)).^2;
            ge_denom = sum(abs(ge_denom).^2);
            temp = max((SNR*ge_num)./(1+ge_denom*SNR));
            ge_temp = [ge_temp (temp)];
        end
        gd = [gd; gd_temp];
        ge = [ge; ge_temp];
    end
    R1 = [];
    for i = 1:14
        R=max(0.5*(log2(1+gd(:,i))-log2(1+ge(:,i))),0);
        R1 = [R1 R];
%         SOP = SOP+Ri;
    end
    R1 = max(R1');
    R1 = R1';
    Ri=(R1<1);
    SOP = SOP+Ri;
end
%%
Pout = SOP/(samples);
semilogy(dB, Pout,'-sr','LineWidth',4);
hold on; grid on;
ylim([10^-4 1]);
xlabel("SNR in dB");
ylabel("SOP");
