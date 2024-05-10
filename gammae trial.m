
%find max value of gamma for each configuration
%for this to happen we need gamma values at every eavesdropper for each configuration
gamma_E = [];
for i = 1:SNR_upperlimit
    gamma_res = [];
    for j = 0:2^M-1
        b = de2bi(j, M);
        denom = 0; num = 0;
        gamma2 = 0;
        for eaves = 1:K
            for temp = 1:M
                if(b(1,temp) == 1) 
                    %temp th node is a jammer
                    denom = denom+G(eaves, temp)*G(eaves, temp)';
                else
                    num = num + G(eaves, temp);
                end
            end
            gamma = SNR*(num*num')/(den+SNR);
            gamma_2 = max(gamma_2, gamma);
        end
        gamma_res = [gamma_res; gamma_2];
    end
    gamma_E = [gamma_E gamma_res];
end
    

%%
SNR = 0;
gamma_temp =[];
for config = 0:2^M-1
b = de2bi(j, M);
for k = 1:K
for i = 1:M
    
    if(b(1, M) == 1)
        den = den + G(i,k)*g(i,k)';
    else 
        num = num + G(i,k);
    end
end
gamma_eavs = (SNR*num*num')/(SNR*den+1);
gamma_max = max(gamma_max, gamma_eavs);
end
gamma_temp = [gamma_temp; gamma_max]
