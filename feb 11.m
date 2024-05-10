clear all; clc;
%% Author: Rohith Kambampati
%% 11 February
h = [0 0.8 0.6 0.7 0 0;
    0 0 0 0.5 0.4 0;
    0 0 0 0.7 0 0.3;
    0 0.5 0.7 0 0.7 0.85;
    0 0.4 0 0.7 0 0.9;
    0 0 0 0 0 0];
e = [0.1 0.15 0.1 0.02;
    0.1 0.3 0.3 0.03;
    0.2 0.15 0.05 0.1;
    0.04 0.05 0.1 0.2;
    0.01 0.04 0.1 0.01;
    0.01 0.03 0.15 0.1];
%% Case 1: Global Defense Cooperative Anti Eavesdropping strategy
A_gd = [];
n_nodes = 6;
n_eaves = 4;
for i = 1 : n_nodes
    for j = 1 : n_nodes
        temp = 0;
        for k = 1 : n_eaves
            temp = temp+e(i,k);
        end
       A_gd(i,j) = temp/h(i,j);
    end
end
%%
for i = 2:5
    if(min(A_gd(i,:))+min(A_gd(:,i)) > 1)
        A_gd(i,:) = Inf;
    end
end
%%
for i = 2:5 
    if(A_gd(i,:) == Inf)
        A_gd(:,i) = Inf;
    end
end
        
%% Case 2: Min-Max Cooperative Anti Eavesdropping strategy
A_mm = [];



