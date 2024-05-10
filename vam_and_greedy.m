clear all; close all; clc;
%%
M = 4; N = 3; K = 2;
var = 1;
e = ones(M,1);
SNR = 10^(0.5);
dB = 0:5:25;
SOP = 0;
samples= 1e6;
th = [0.5 0.6 0.7 0.8 0.9];
for iii = 1:samples
    gd = []; ge = [];
    H = sqrt(var/2)*(randn(N,M)+i*randn(N,M));
    G = sqrt(var/2)*(randn(K,M)+i*randn(K,M));
    all_nodes = [1 1 1 1]; % start with assumption as all of them as relays and w
    Q = H*e*e'*H';
    [u, v] = eig(Q);
    [~, label] = max(max(abs(v)));
    w = u(:,2);
    [~, max_label] = max(abs(w' * H));
    c = [];
    for i = 1:M
        temp = abs(H(:,max_label)'*H(:,i))/(abs(H(:,max_label)')*abs(H(:,i)));
        c = [c temp];
    end
    [~, min_c] = min(c);
    relay_set = (c>0)*1;
    jammer_set = 1-relay_set;
    if(jammer_set == [0 0 0 0])
        [~,relay_min] = min(relay_set);
        jammer_set(relay_min) = 1;
        relay_set(relay_min) = 0;
    end
    
    H_J = H*diag(jammer_set);
    H_R = H*diag(relay_set);
    
    for ii = 1:length(dB)
        SNR = 10^(dB(ii)/10);
        Q = inv(SNR*H_J*H_J' + eye(N)) * SNR*(H_R*e)*(H_R*e)';
        gd_temp = max(abs(eig(Q)));
        
        [~,L2] = find(jammer_set==1);
        [~,L1] = find(relay_set==1);
        ge_num = G(:, L1);
        ge_denom = G(:, L2);
        ge_num = abs(sum(ge_num)).^2;
        ge_denom = sum(abs(ge_denom).^2);
        ge_temp = max(SNR*ge_num./(1+ge_denom*SNR));
        
        gd = [gd gd_temp];
        ge = [ge ge_temp];
    end
    R = max(0.5* (log2(1+gd)-log2(1+ge)),0);
    Ri = (R<1);
    SOP = SOP+Ri;
end
%%
Pout = SOP/samples;
semilogy(dB, Pout,'-sr','LineWidth',2);
hold on; grid on;
ylim([10^-6 1]);
xlabel("SNR in dB");
ylabel("SOP");
hold on;
%%
clear all; clc;
%% Generating H and G matrices
M = 4; N = 3; K = 2;
var = 1;
e = ones(M,1);
dB = 0:5:25;
samples= 1e4;
SOP=0;
for iii = 1:samples
    gd1 = []; ge1 = [];
    gd2 = []; ge2 = [];
    gd3 = []; ge3 = [];
    H = sqrt(var/2)*(randn(N,M)+i*randn(N,M));
    G = sqrt(var/2)*(randn(K,M)+i*randn(K,M));
    for ii = 1:length(dB)
        SNR = 10^(dB(ii)/10);
        relay_set = ones(M,1);
        % case : L2 = 0
        Q = (H*e)*(H*e)';
        [u1, v1] = eig(Q);
        [~, label] = max(abs(eig(Q)));
        w = u1(:, label); % finding max eigen vector
        temp = w' * H;
        j = min(abs(temp));
        [~, jammer1] = find(abs(temp) == j);% finding node with lease w' * h
        jammer_set = zeros(M,1);
        jammer_set(jammer1) = 1;
        relay_set = 1-jammer_set;% now L2 is made 1 and L1 = 3
        
        H_R = H*diag(relay_set);
        H_J = H*diag(jammer_set);
        Q = inv(SNR*H_J*H_J' + eye(N)) * SNR*(H_R*e)*(H_R*e)';
        gd_temp = max(abs(eig(Q)));
        
        [~,L2] = find(jammer_set==1);
        [~,L1] = find(relay_set==1);
        ge_num = G(:, L1);
        ge_denom = G(:, L2);
        ge_num = abs(sum(ge_num)).^2;
        ge_denom = sum(abs(ge_denom).^2);
        ge_temp = max(SNR*ge_num./(1+ge_denom*SNR));
        
        gd1 = [gd1 gd_temp];
        ge1 = [ge1 ge_temp];
        %%
        % case : L2 = 1
        b = diag(jammer_set);
        H_J = H*b;
        H_R = H*(1-b);
        Q = inv(SNR*H_J*H_J' + eye(N)) * SNR*(H_R*e)*(H_R*e)';
        [u2, v2] = eig(Q);
        [~, label] = max(abs(eig(Q)));
        w = u2(:, label); % finding max eigen vector
        temp = w'*H_R;
        temp1 = abs(temp);
        temp1(jammer1) = 1000;
        j = min(abs(temp1));
        
        [~, jammer2] = find(abs(w' * H_R) == j);% finding node with lease w' * h
        jammer_set(jammer2(1,length(jammer2))) = 1;% If there are more than one node with min w'* H, then we choose the one closest to the destination as the jammer
        relay_set = 1-jammer_set;
        
        H_R = H*diag(relay_set);
        H_J = H*diag(jammer_set);
        Q = inv(SNR*H_J*H_J' + eye(N)) * SNR*(H_R*e)*(H_R*e)';
        gd_temp = max(abs(eig(Q)));
        
        [~,L2] = find(jammer_set==1);
        [~,L1] = find(relay_set==1);
        ge_num = G(:, L1);
        ge_denom = G(:, L2);
        ge_num = abs(sum(ge_num)).^2;
        ge_denom = sum(abs(ge_denom).^2);
        ge_temp = max(SNR*ge_num./(1+ge_denom*SNR));
        
        gd2 = [gd2 gd_temp];
        ge2 = [ge2 ge_temp];
        %%
%       case : L2 = 2
        b = diag(jammer_set);
        H_J = H*b;
        H_R = H*(1-b);
        Q = inv(SNR*H_J*H_J' + eye(N)) * SNR*(H_R*e)*(H_R*e)';
        [u3, v3] = eig(Q);
        [~, label] = max(abs(eig(Q)));
        w = u3(:, label); % finding max eigen vector
        temp = w'*H_R;
        temp1 = abs(temp);
        for i = 1:length(jammer_set)
            if(jammer_set(i, 1) == 1)
                temp1(1, i) = 1000;
            end
        end
        [~, jammer3] = min(temp1);
        jammer_set(jammer3(1,length(jammer3))) = 1;
        relay_set = 1 - jammer_set;
        
        %%
        H_R = H*diag(relay_set);
        H_J = H*diag(jammer_set);
        Q = inv(SNR*H_J*H_J' + eye(N)) * SNR*(H_R*e)*(H_R*e)';
        gd_temp = max(abs(eig(Q)));
        
        [~,L2] = find(jammer_set==1);
        [~,L1] = find(relay_set==1);
        ge_num = G(:, L1);
        ge_denom = G(:, L2);
        ge_num = abs(sum(ge_num)).^2;
        ge_denom = sum(abs(ge_denom).^2);
        ge_temp = max(SNR*ge_num./(1+ge_denom*SNR));
        
        gd3 = [gd3 gd_temp];
        ge3 = [ge3 ge_temp];
    end
    R1=max(0.5* (log2(1+gd1)-log2(1+ge1)),0);
    R2=max(0.5* (log2(1+gd2)-log2(1+ge2)),0);
    R3=max(0.5* (log2(1+gd3)-log2(1+ge3)),0);
    Ri=(max(max(R1, R2), R3) <1 );
    SOP = SOP+Ri;
end
%%
Pout = SOP/samples;
semilogy(dB, Pout,'-sr','LineWidth',2);
hold on; grid on;
ylim([10^-6 1]);
xlabel("SNR in dB");
ylabel("SOP");
legend("VAM", "Greedy");
%%


