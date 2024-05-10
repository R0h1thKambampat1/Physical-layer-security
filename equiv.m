%%%%%%%%%%%%%5
%%% code for paper: Zhenguo Gao, Danjie Chen, etc. Outage probability equivalency of three typical relay 
%%% selection schemes for selective DF relay networks with selection combining. 
%%% Wireless Personal Communications (Springer). Jun, 2015, 85(3):1205-1215 
%%%%%%%%%%%%%%
function []=equiv()
clear all;
format long;
rand('state',sum(100*clock));
digits=64;

scn=10^6;
num_sim=1;
num_testscn=1;
        
P_max=0.5;
hp_min=1;
hp_max=10;

global rth_db;
rth_db=0;
rth=10^(rth_db/10);
N0=1;

list_snrdb=-10:5:20;

ii_ProRS=1;
ii_ReRS=2;
ii_BReRS=3;

num_p=1;
num_r=3;

su=[0,0];
du=[1,0];
pu(:,1)=0.15+0.15*[1:1:num_p]';
pu(:,2)=0.3;
ru=zeros(num_r,2);
ru(:,1)=0.4;
ru(:,2)=0.1+1/num_r.*([0:1:num_r-1]'-floor(num_r/2));

if 1,
   figure; %show network topology;
   hold on;
   plot(su(1),su(2),'ro');
   plot(du(1),du(2),'rx');
   plot(pu(:,1)',pu(:,2)','g*');
   plot(ru(:,1)',ru(:,2)','bd');
   axis([-0.1,1.1,-1,1]);
   box on;
   grid on;
   title('Network Topology');
end

lsd=sqrt((du(1)-su(1))^2+(du(2)-su(2))^2);
lsp=sqrt((pu(:,1)-su(1)).^2+(pu(:,2)-su(2)).^2);
lsr=sqrt((ru(:,1)-su(1)).^2+(ru(:,2)-su(2)).^2);
lrd=sqrt((ru(:,1)-du(1)).^2+(ru(:,2)-du(2)).^2);

for i=1:1:num_r,
   lrp(:,i)=sqrt((ru(i,1)*ones(num_p,1)-pu(:,1)).^2+(ru(i,2)*ones(num_p,1)-pu(:,2)).^2);
end

alpha=-4;
mat_gsd=lsd^alpha;
mat_gsr=lsr.^alpha;
mat_grd=lrd.^alpha;
mat_gsp=lsp.^alpha;
mat_grp=lrp.^alpha;

list_snr=10.^(list_snrdb/10);
snr_num=length(list_snr);

mat_out=zeros(3,snr_num);

mat_IDr(1,1:1)=[3];
mat_IDr(2,1:1)=[3];
mat_IDr(3,1:3)=[1,2,3];
mat_IDr(4,1:3)=[1,2,3];

gsd=mat_gsd;
gsr=mat_gsr;
grd=mat_grd;
gsp=mat_gsp;
grp=mat_grp;

list_opt_all=[];
list_opt_sr=[];
list_opt_rd=[];

cdf_ploted=0;

for ii_test=1:1:num_testscn,
   
   list_IDr=mat_IDr(ii_test,find(mat_IDr(ii_test,:)>0.5));
   num_IDr=length(list_IDr);
    
   for ii_snr=1:1:snr_num,
      Ip=list_snr(ii_snr);
      Psdb=list_snrdb(ii_snr);
      for ii_sim=1:1:num_sim,
         
         % progress indication
         step_str=strcat('[Scn(',num2str(num_testscn),':',num2str(ii_test),')]');
         step_str=strcat(step_str,'[SNR(',num2str(snr_num),':',num2str(ii_snr),')]');
         step_str=strcat(step_str,'[Sim(',num2str(num_sim),':',num2str(ii_sim),')]');
         disp(step_str);
         
         hsd2=exprnd(gsd,1,scn);
         hsr2=exprnd(gsr*ones(1,scn));
         hsp2=exprnd(gsp,1,scn);
         hrp2=exprnd(grp'*ones(1,scn));

         
         hsp2=hp_min+(hp_max-hp_min)*rand(1,scn);
         hrp2=hp_min+(hp_max-hp_min)*rand(num_r,scn);

         for i=1:1:num_r,
            hrd2(i,:)=grd(i)*rand_rician(0,10,scn);
         end
         
         Ps=Ip./hsp2;
         Pr=Ip./hrp2;
         
         %Ps=min(P_max,Ps);
         %Pr=min(P_max,Pr);
         
         snr_sd_mat=Ps.*hsd2/N0;
         snr_sr_mat=ones(num_r,1)*Ps.*hsr2/N0;
         snr_rd_mat=Pr.*hrd2/N0;
         
         snr_sd_mat_temp=zeros(1,scn);
         
         % all optimal
         snr_min_r_mat=[];
         for ii_r=1:1:num_IDr,      
            id_r=list_IDr(ii_r);
            snr_min_r_mat(i,:)=min([snr_sr_mat(id_r,:);snr_rd_mat(id_r,:)],[],1);
         end
         snr_opt_all=max(snr_min_r_mat,[],1);
         
         % rd optimal
         snr_srsuc_rdmax=zeros(num_IDr,scn);
         for ii_r=1:1:num_IDr,      
            id_r=list_IDr(ii_r);
            suc_sr_list=find((snr_sr_mat(id_r,:))>=rth);
            snr_srsuc_rdmax(i,suc_sr_list)=snr_rd_mat(id_r,suc_sr_list);
         end
         snr_opt_rd=max(snr_srsuc_rdmax,[],1);
         
         % sr optimal
         snr_rdsuc_srmax=zeros(num_IDr,scn);
         for ii_r=1:1:num_IDr,      
            id_r=list_IDr(ii_r);
            suc_rd_list=find((snr_rd_mat(id_r,:))>=rth);
            snr_rdsuc_srmax(i,suc_rd_list)=snr_sr_mat(id_r,suc_rd_list);
         end
         
         snr_opt_sr=max(snr_rdsuc_srmax,[],1);
         
         snr_opt_all=max([snr_opt_all;snr_sd_mat],[],1);
         snr_opt_sr=max([snr_opt_sr;snr_sd_mat],[],1);
         snr_opt_rd=max([snr_opt_rd;snr_sd_mat],[],1);
         
         if cdf_ploted==0,
            [outnum]=length(find(snr_opt_all<rth));
            snr_opt=[snr_opt_all;snr_opt_sr;snr_opt_rd];
            calAndPlotCDF(snr_opt,Psdb,outnum);
            %cdf_ploted=1;
            %pause;
         end
         
         
         mat_out(ii_ProRS,ii_snr)=mat_out(ii_ProRS,ii_snr)+length(find(snr_opt_all<rth));
         mat_out(ii_ReRS,ii_snr)=mat_out(ii_ReRS,ii_snr)+length(find(snr_opt_sr<rth));
         mat_out(ii_BReRS,ii_snr)=mat_out(ii_BReRS,ii_snr)+length(find(snr_opt_rd<rth));
         
      end %for sim
   end %for snr
end %for test

mat_out=mat_out/num_sim/scn;
mat_out;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;


semilogy(list_snrdb, mat_out(ii_ProRS,:),'rs','MarkerEdgeColor','r', 'MarkerFaceColor','w', 'MarkerSize',10);
hold on;
semilogy(list_snrdb, mat_out(ii_ReRS,:),'b*');
semilogy(list_snrdb, mat_out(ii_BReRS,:),'-k');

str_legend={'ProRS','ReRS','BReRS'};  

mat_out

legend(str_legend);
xlabel('Interference power limit {\it{I_P}}(dB)');
ylabel('Outage probability');
set(gca,'YMinorGrid','on');
set(gca,'XGrid','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555
%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55%5

function [rnd_v]=rand_rayleign(type_rand,num)
type_power=0;
type_gain=1;
type_signal=2;

type_rand=type_power;

x1=randn(1,num)/sqrt(2);
x2=randn(1,num)/sqrt(2);
rnd_v=zeros(1,num);
if type_rand==type_power
   rnd_v=x1.^2+x2.^2;
elseif type_rand==type_gain
   rnd_v=sqrt(x1.^2+x2.^2);
elseif type_rand==type_signal
   rnd_v=x1+j.*x2;
else
   disp('unknow rand type!!!');
   return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555
%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55%5
function [rnd_v]=rand_rician(type_rand, k,num)
type_power=0;
type_gain=1;
type_signal=2;

type_rand=type_power;

x1=randn(1,num)/sqrt(2);
x2=randn(1,num)/sqrt(2);

rnd_v=zeros(1,num);
sqrtk=sqrt(k);

if type_rand==type_power
   rnd_v=((sqrtk+x1).^2+x2.^2)/(k+1);
elseif type_rand==type_gain
   rnd_v=sqrt((sqrtk+x1).^2+x2.^2)/sqrt(k+1);
elseif type_rand==type_signal
   rnd_v=(sqrtk+x1)+j.*x2/sqrt(k+1);
else
   disp('unknow rand type!!!');
   return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555
%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55%5
function []=calAndPlotCDF(snr_opt,Psdb,outnum)
[curve_num,snr_num]=size(snr_opt);
[snr_value,snr_ii]=sort(snr_opt,2);
figure;
x_list=(1:1:snr_num)/snr_num*100;

colorstr={'r','b','k','g','m','c','y','r','g','b','k','y'};
markstylestr={'s','d','v','^','>','<','p','h','o','x','+','*',};
linestylestr={'-',':','-.','--','-',':','-.','--','-',':','-.','--'};
errorlinestylestr={'-r','-b','-k','-g','-m','-c','-y','-r','-g','-b','-k','-y'};
global rth_db;

for i=1:1:curve_num,
   plot(snr_value(i,:),x_list,strcat(char(linestylestr(i)),char(colorstr(i))));
   hold on;
end

for i=1:1:curve_num,
   if outnum==0,
       disp('outnum=0');
   end
   plot(snr_opt(i,snr_ii(i,outnum)),outnum/snr_num*100,strcat(char(linestylestr(i)),char(colorstr(i)),char(markstylestr(i))), 'lineWidth', 1,'MarkerEdgeColor',char(colorstr(i)),'MarkerFaceColor','w', 'MarkerSize',5);
end

axis([0  10  -Inf Inf]); 

str_legend={'ProRS','ReRS','BReRS'};  
legend(str_legend);
xlabel('end-to-end SNR \gamma_{E2E}');
ylabel('F_{\gamma_{E2E}}(\gamma) (%)');
title(strcat('F_{\gamma_{E2E}}(\gamma)(\gamma_{th}=',num2str(rth_db),'dB, N_0=1,I_P=',num2str(Psdb),'dB)'));

