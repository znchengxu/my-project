

function [queding_houxuan,yuanxian_houxuan,Reconstruction_rate,APRESS]=modeling_term(y,houxuan)
    for i=1:length(houxuan)
        data=[];
        data=houxuan{i,2};
        xi_m(i)=(y*data').^2/((y*y')*(data*data'));
    end

    [shu_n,shu_m]=max(xi_m);
    queding_houxuan=houxuan(shu_m,:);
    yuanxian_houxuan=houxuan(shu_m,:);

    shuju=[];
    for xs=1:size(queding_houxuan,1)
        shuju(xs,:)=yuanxian_houxuan{xs,2};
    end
    canshu=[];
    canshu=pinv(shuju*(shuju)')*shuju*y';
    APRESS(length(canshu))=mean((shuju*canshu-y).^2);
    Reconstruction_rate=1-(var(y-canshu*shuju))/var(canshu*shuju);
    houxuan(shu_m,:)=[];

    while size(queding_houxuan,1)<=40
        num=size(queding_houxuan,1)+1;
        zhengjiao_houxuan=[];
        zhengjiao_houxuan=houxuan;        
        for i=1:length(houxuan)   
            pp=zeros(size(y,1),size(y,2)); 
            for k=1:num-1                                     
                pp=pp+(houxuan{i,2}*queding_houxuan{k,2}')/(queding_houxuan{k,2}*queding_houxuan{k,2}').*(queding_houxuan{k,2});
            end        
            zhengjiao=[];           
            zhengjiao=houxuan{i,2}-pp;             
            zhengjiao_houxuan{i,2}=zhengjiao;         
        end
        xi_m=[];
        for i=1:length(zhengjiao_houxuan)
            data=[];
            data=zhengjiao_houxuan{i,2};
            xi_m(i)=(y*data').^2/((y*y')*(data*data'));
        end
        [shu_n,shu_m]=max(xi_m);
        queding_houxuan(num,:)=zhengjiao_houxuan(shu_m,:);
        yuanxian_houxuan(num,:)= houxuan(shu_m,:);
        houxuan(shu_m,:)=[];
        shuju=[];
        for xs=1:length(queding_houxuan)
            shuju(xs,:)=yuanxian_houxuan{xs,2};
        end
        canshu=[];
        canshu=pinv(shuju*(shuju)')*shuju*y';
%         APRESS(length(canshu))=(1/(1-(length(canshu))/(length(y)))).^2*1/length(y)*mean((shuju'*canshu-y').^2);
        APRESS(length(canshu))=mean((canshu'*shuju-y).^2);
        Reconstruction_rate(num)=1-(var(y-canshu'*shuju))/var(canshu'*shuju);
    end
end