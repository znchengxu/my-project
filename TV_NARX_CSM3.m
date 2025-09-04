clear
clc
tic
% load data
load('C:\Users\lab\Desktop\稀疏建模\Benchmark_EEG_small\Benchmark_EEG_small.mat')
U=data.u;
Y=data.y;  

%model structure detection for each subject
for ren=1:size(U,1) %
    y=[];
    u=[];
    for i=1:size(U,2)-1%modeling by the first 6 sets of data
        u1=[];
        y1=[];
        u1=squeeze(U(ren,i,:));
        y1=squeeze(Y(ren,i,:)); 
        u1=mapminmax(u1',-1,1);
        y1=mapminmax(y1',-1,1);
        if i==1
            u=u1;
            y=y1;
        else
            u=[u u1];
            y=[y y1];
        end
    end

    houxuan=[];
    houxuan=dictionaries_creater(y,u,5,20,2);
    final_houxuan(ren)={houxuan};

    queding_houxuan=[];yuanxian_houxuan=[];Reconstruction_rate=[];
    [queding_houxuan,yuanxian_houxuan,Reconstruction_rate,APRESS]=modeling_term(y,houxuan);
    final_queding_houxuan(ren)={queding_houxuan};
    final_yuanxian_houxuan(ren)={yuanxian_houxuan};
    final_Reconstruction_rate(ren)={Reconstruction_rate};
    final_APRESS(ren)=APRESS(end);
    final_y(ren)={y};
end

%% Distinguish between the shared terms and non-shared terms
num_chong=1;
for i=1:size(final_yuanxian_houxuan{1, 1},1)
    xiang=[];
    xiang=final_yuanxian_houxuan{1, 1}{i, 1};
    num=0;
    for j=2:length(final_yuanxian_houxuan)
        data=[];
        data=final_yuanxian_houxuan{1,j};
        for k=1:length(data)
            if strcmp(xiang,data{k,1})
                num=num+1;
                break;
            end
        end
    end
    if num==length(final_yuanxian_houxuan)-1 
        chongfu_term(num_chong)={xiang};num_chong=num_chong+1;%the shared terms
    end
end


term=houxuan(:,1);
weizhi_chongfu=[];
for j=1:length(chongfu_term)
    ax=[];
    ax=chongfu_term{j};
    for k=1:size(houxuan,1)
        if strcmp(ax,term{k,1})
            weizhi_chongfu(j)=k;
        end
    end        
end

%Count the number of each non-shared regression term
for i=1:length(final_yuanxian_houxuan)
    data=[];
    data=final_yuanxian_houxuan{1,i};
    num=1;
    for j=1:length(data)
        if ismember(data{j,1}, chongfu_term)
            jilu(num)=j;
            num=num+1;
        end
    end 
    data(jilu,:)=[];
    feichongfu{1,i}={};
    feichongfu{1,i}=data(:,1);
end

for i=1:length(feichongfu)
    data=[];
    data=feichongfu{1,i};
    if i==1
        feichongfu_term1=data;
    else
        feichongfu_term1(length(feichongfu_term1)+1:length(feichongfu_term1)+length(data))=data;
    end
end
feichongfu_term=unique(feichongfu_term1);%the non-shared terms

cishu_feichongfu_term=zeros(1,length(feichongfu_term));
for i=1:length(feichongfu_term)
    data=[];
    data=feichongfu_term{i,1};
    num=0;
    for j=1:length(feichongfu)
        data1=[];
        data1=feichongfu{1,j};
        for k=1:length(data1)
            if strcmp(data1{k,1}, data)
                num=num+1;
                break
            end
        end
    end
    cishu_feichongfu_term(i)=num;
end
yuzhi=1-cishu_feichongfu_term./size(U,1);%Set the selection threshold for non-shared items


%% Use genetic algorithms to select non-shared terms

geshu=25;
NIND=100;
MAXGEN=100;
GGAP=0.90;
PX=0.8;
PM1=0.01;
PM2=PM1;PM3=PM1;
trace=zeros(MAXGEN,length(yuzhi)+1);
zhongqun=zeros(NIND,length(feichongfu_term));
for i=1:size(zhongqun,1)
%     if i<NIND/2
%         k=0;
%         while k<1
%             bu=[];
%             bu=randperm(size(zhongqun,2),1);
%             if rand(1)>yuzhi(bu)
%                 zhongqun(i,bu)=1;
%             end
%             if length(find(zhongqun(i,:)==1))==15%指向性
%                 k=1;
%             end
%         end    
%     else
        bu=[];           
        bu=randperm(size(zhongqun,2),geshu);                    
        zhongqun(i,bu)=1;             
%     end
end

%Calculate the fitness
zong_fitness=fit(zhongqun,feichongfu_term,houxuan,weizhi_chongfu,final_houxuan,final_y);

%The populations are divided into three groups based on fitness.
[zong_fitness_n,zong_fitness_mm]=sort(zong_fitness);
zhongqun_h=zhongqun(zong_fitness_mm(1:34),:);
zong_fitness_h=zong_fitness(zong_fitness_mm(1:34));
zhongqun_m=zhongqun(zong_fitness_mm(35:34+33),:);
zong_fitness_m=zong_fitness(zong_fitness_mm(35:67));
zhongqun_s=zhongqun(zong_fitness_mm(68:end),:);
zong_fitness_s=zong_fitness(zong_fitness_mm(68:end));


gen=0;
while gen<MAXGEN
    %
    zong_fitness_n=[];zong_fitness_mm=[];
    [zong_fitness_n,zong_fitness_mm]=sort(zong_fitness);
    zhongqun_h=zhongqun(zong_fitness_mm(1:34),:);
    zong_fitness_h=zong_fitness(zong_fitness_mm(1:34));
    zhongqun_m=zhongqun(zong_fitness_mm(35:34+33),:);
    zong_fitness_m=zong_fitness(zong_fitness_mm(35:67));
    zhongqun_s=zhongqun(zong_fitness_mm(68:end),:);
    zong_fitness_s=zong_fitness(zong_fitness_mm(68:end));
    chushi_var1=var(zong_fitness_h);
    chushi_var2=var(zong_fitness_m);
    chushi_var3=var(zong_fitness_s);
% chushi_var4=var(zong_fitness_hunhe);
%     suiji=[];
%     suiji=randperm(33,11);
%     zong_fitness_hunhe=zong_fitness_h(suiji);
%     zong_fitness_hunhe(12:11*2)=zong_fitness_m(suiji);
%     zong_fitness_hunhe(23:33)=zong_fitness_s(suiji);
%     zhongqun_hunhe=zhongqun_h(suiji,:);
%     zhongqun_hunhe(12:11*2,:)=zhongqun_m(suiji,:);
%     zhongqun_hunhe(23:33,:)=zhongqun_s(suiji,:);
    %

   %Selection, crossover, mutation
   FitnV1=ranking(zong_fitness_h');
   SelCh1=select('sus',zhongqun_h,FitnV1,GGAP); 
   SelCh1=recombin('xovsp',SelCh1,PX);
   SelCh1=mut(SelCh1,PM1);
   SelCh1=tran(SelCh1,geshu,zong_fitness_h,zhongqun_h);
    
   FitnV2=ranking(zong_fitness_m');
   SelCh2=select('sus',zhongqun_m,FitnV2,GGAP); 
   SelCh2=recombin('xovsp',SelCh2,PX);
   SelCh2=mut(SelCh2,PM2);
   SelCh2=tran(SelCh2,geshu,zong_fitness_m,zhongqun_m);
    
   FitnV3=ranking(zong_fitness_s');
   SelCh3=select('sus',zhongqun_s,FitnV3,GGAP); 
   SelCh3=recombin('xovsp',SelCh3,PX);
   SelCh3=mut(SelCh3,PM3);
   SelCh3=tran(SelCh3,geshu,zong_fitness_s,zhongqun_s);
  

   zong_fitness11=[];
   zong_fitness11=fit(SelCh1,feichongfu_term,houxuan,weizhi_chongfu,final_houxuan,final_y);
   zong_fitness12=[];
   zong_fitness12=fit(SelCh2,feichongfu_term,houxuan,weizhi_chongfu,final_houxuan,final_y);
   zong_fitness13=[];
   zong_fitness13=fit(SelCh3,feichongfu_term,houxuan,weizhi_chongfu,final_houxuan,final_y);

   
   %Population update 
   Z_F=[];
   Z_F=[zong_fitness_h zong_fitness11];
   Z_F_zhongqun=[];
   Z_F_zhongqun=[zhongqun_h ;SelCh1];
   Z_F_n=[];Z_F_m=[];
   [Z_F_n,Z_F_m]=sort(Z_F);
   zong_fitness_h=[];
   zong_fitness_h=Z_F(Z_F_m(1:size(zhongqun_h,1)));
   zhongqun_h=[];
   zhongqun_h=Z_F_zhongqun(Z_F_m(1:size(zong_fitness_h,2)),:);
   %
   Z_F=[];
   Z_F=[zong_fitness_m zong_fitness12];
   Z_F_zhongqun=[];
   Z_F_zhongqun=[zhongqun_m ;SelCh2];
   Z_F_n=[];Z_F_m=[];
   [Z_F_n,Z_F_m]=sort(Z_F);
   zong_fitness_m=[];
   zong_fitness_m=Z_F(Z_F_m(1:size(zhongqun_m,1)));
   zhongqun_m=[];
   zhongqun_m=Z_F_zhongqun(Z_F_m(1:size(zong_fitness_m,2)),:);
   % 
   Z_F=[];
   Z_F=[zong_fitness_s zong_fitness13];
   Z_F_zhongqun=[];
   Z_F_zhongqun=[zhongqun_s ;SelCh3];
   Z_F_n=[];Z_F_m=[];
   [Z_F_n,Z_F_m]=sort(Z_F);
   zong_fitness_s=[];
   zong_fitness_s=Z_F(Z_F_m(1:size(zhongqun_s,1)));
   zhongqun_s=[];
   zhongqun_s=Z_F_zhongqun(Z_F_m(1:size(zong_fitness_s,2)),:);
    % 
   
   
%    suiji=[];
%    suiji=randperm(33,2);
%    zancun1=[];
%    zancun2=[];
%    zancun3=[];
% %    zancun4=[];
%    zc1=[];
%    zc2=[];
%    zc3=[];
% %    zc4=[];
%    zancun1=zhongqun_h(suiji,:);
%    zancun2=zhongqun_m(suiji,:);
%    zancun3=zhongqun_s(suiji,:);
% %    zancun4=zhongqun_hunhe(suiji,:);
%    zc1=zong_fitness_h(suiji);
%    zc2=zong_fitness_m(suiji);
%    zc3=zong_fitness_s(suiji);
% %    zc4=zong_fitness_hunhe(suiji);
%    zhongqun_h(suiji(1),:)=zancun2(1,:);
%    zhongqun_h(suiji(2),:)=zancun3(1,:);
% %    zhongqun_h(suiji(3),:)=zancun4(1,:);
%    
%    zhongqun_m(suiji(1),:)=zancun1(1,:);
%    zhongqun_m(suiji(2),:)=zancun3(2,:);
% %    zhongqun_m(suiji(3),:)=zancun4(2,:);
%    
%    zhongqun_s(suiji(1),:)=zancun1(2,:);
%    zhongqun_s(suiji(2),:)=zancun2(2,:);
% %    zhongqun_s(suiji(3),:)=zancun4(3,:); 
%    
% %    zhongqun_hunhe(suiji(1),:)=zancun1(3,:);
% %    zhongqun_hunhe(suiji(2),:)=zancun2(3,:);
% %    zhongqun_hunhe(suiji(3),:)=zancun3(3,:);
%    
%    zong_fitness_h(suiji(1))=zc2(1);
%    zong_fitness_h(suiji(2))=zc3(1);
% %    zong_fitness_h(suiji(3))=zc4(1);
%    
%    zong_fitness_m(suiji(1))=zc1(1);
%    zong_fitness_m(suiji(2))=zc3(2);
% %    zong_fitness_m(suiji(3))=zc4(2);
%    
%    zong_fitness_s(suiji(1))=zc1(2);
%    zong_fitness_s(suiji(2))=zc2(2);
% %    zong_fitness_s(suiji(3))=zc4(3);
%    
% %    zong_fitness_hunhe(suiji(1))=zc1(3);
% %    zong_fitness_hunhe(suiji(2))=zc2(3);
% %    zong_fitness_hunhe(suiji(3))=zc3(3);
   
   %Adjust the coefficient of mutation
   muqian_var1=var(zong_fitness_h);
   muqian_var2=var(zong_fitness_m);
   muqian_var3=var(zong_fitness_s);
   if muqian_var1<chushi_var1*0.1
       PM1=PM1*1.05;
   end
    if muqian_var2<chushi_var2*0.1
       PM2=PM1*1.05;
    end
    if muqian_var3<chushi_var3*0.1
       PM3=PM1*1.05;
   end
    
 
   gen=gen+1;
   a=[];b=[];
   zong_fitness1=[];
   zong_fitness1=[zong_fitness_h,zong_fitness_m,zong_fitness_s];
   zong_fitness1_n=[];zong_fitness1_m=[];
   [zong_fitness1_n,zong_fitness1_m]=sort(zong_fitness1);
   zong_fitness=[];
   zong_fitness=zong_fitness1_n(1:100);
   %
   zhongqun1=[];
   zhongqun1=[zhongqun_h;zhongqun_m;zhongqun_s];
   zhongqun=[];
   zhongqun=zhongqun1(zong_fitness1_m(1:100),:);
   [a,b]=min(zong_fitness);
   trace(gen,1:size(zhongqun,2))=zhongqun(b,:);
   trace(gen,size(zhongqun,2)+1)=a;
end

%Save the optimal structure
term=houxuan(:,1);
xz=trace(gen,1:size(trace,2)-1);
ax=[];
ax=feichongfu_term(find(xz==1));
for i=1:length(ax)
    for k=1:size(houxuan,1)
        if strcmp(ax{i,1},term{k,1})
            weizhi(i)=k;
        end 
    end
end
ZONG_weizhi=[weizhi_chongfu weizhi];
for ren=1:10
    dy=[];
    dy=final_houxuan{1,ren};
    xuanze{1,ren}=dy(ZONG_weizhi,:);
end

%Solve the parameters of the optimal sharing structure
for i=1:length(xuanze)
    zancun={};
    zancun=xuanze{1,i}; 
    shuju=[];
    for k=1:length(zancun)
        shuju(k,:)=zancun{k,2};
    end
    mubiao_y=[];
    mubiao_y=final_y{1, i};
    canshu=[];
    canshu=pinv(shuju*(shuju)')*shuju*mubiao_y';    
    zuizhongcanshu{1,i}=canshu;
    akk(i)=mean((canshu'*shuju-mubiao_y).^2);
end
term=houxuan(:,1);
num=1;
for i=1:length(xuanze{1,1})
    baoliuxiang=[];
    baoliuxiang=xuanze{1, 1}{i,1};
    for j=1:length(term)
         if strcmp(baoliuxiang,term{j,1})
             jilu_weizhi1(num)=j;
             num=num+1;
         end
    end
end

%% prediction
for zu=1:10
    
    y=[];u=[];
    u=squeeze(U(zu,7,:));
    y=squeeze(Y(zu,7,:));
    u=mapminmax(u',-1,1);
    y=mapminmax(y',-1,1);
    
    houxuan=[];
    houxuan=dictionaries_creater(y,u,5,20,2);
    hxx=[];
    for i=1:length(jilu_weizhi1)
        hxx(i,:)=houxuan{jilu_weizhi1(i),2};
    end
    TH=[];TH=zuizhongcanshu{1,zu};
    Pre_y1=TH'*hxx;
%     Pre_y1(1:20)=y(1:20);
%     CC1=cov(y,Pre_y1)/(std(y)*std(Pre_y1));
%     CC=CC1(1,2);
%     VAF=[];
%     VAF=1-(var(y-Pre_y1))/var(Pre_y1);
%     jieguo(zu)=VAF;
    MSE(zu)=mean((y-Pre_y1).^2);
end
sum(akk)
sum(MSE)
toc