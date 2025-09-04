%% Dictionary generation
 function houxuan=dictionaries_creater(y,u,ny,nx,order)

n=length(y);

%The lag term of y
for i=1:ny
    dairu=[];
    dairu=zeros(1,i);
    YY=y;
    YY=[dairu YY];
    YY(n+1:end)=[];
    a=[];
    a = ['y', num2str(i)];
    houxuan_y(i,1)={a};
    houxuan_y(i,2)={YY};
end

%The lag term of x
if nx==0
    houxuan_p=[];
else
    for i=1:nx
        dairu=[];
        dairu=zeros(1,i);
        XX=u;
        XX=[dairu XX];
        XX(n+1:end)=[];
        a=[];
        a = ['p', num2str(i)];
        houxuan_p(i,1)={a};
        houxuan_p(i,2)={XX};
    end
end

huizong=[houxuan_p;houxuan_y];
A=size(huizong,1);
A_str = num2str(A); 
num_digits = numel(A_str); 

%Obtain candidate terms
houxuan=[];
for i=1:length(huizong(:,1))
    zancun=[];jieshu=1;
    zancun=huizong(i,:);
    if numel(houxuan)==0
         houxuan=zancun;
    else
        houxuan(length(houxuan(:,1))+1,:)=zancun;
    end
    zhongjian_houxuan=[];
    zhongjian_houxuan=zancun;
    while jieshu<order% 
        zan_houxuan=[];
        for k=1:length(zhongjian_houxuan(:,1))
            su=[];
            su=zhongjian_houxuan(k,:);
            for j=1:length(huizong)
                zan_houxuan{j,1}=[su{1,1},huizong{j,1}];
                zan_houxuan(j,2)={su{1,2}.*huizong{j,2}};
            end
            houxuan(length(houxuan(:,1))+1:length(houxuan(:,1))+length(zan_houxuan(:,1)),:)=zan_houxuan;          
        end
        zhongjian_houxuan=[];
        zhongjian_houxuan=zan_houxuan;
        jieshu=jieshu+1;
    end
end

%Eliminate duplicate candidate terms
mingcheng=houxuan(:,1);
for i=1:length(mingcheng)
    str=mingcheng{i,1};  
    letters=[];
    letters = regexp(str, '\D', 'match'); 
    weizhi=[];
    for j=1:length(letters)
        weizhi1=strfind(str,letters{1,j});
        weizhi=[weizhi weizhi1];
        weizhi=unique(weizhi);
    end
    cucu=[];
    for j=1:length(weizhi)
        if j==length(weizhi)
            mingcheng1={str(weizhi(j):length(str))};
        else       
            mingcheng1={str(weizhi(j):weizhi(j+1)-1)};
        end
        WEIZHI=ismember(huizong(:,1),mingcheng1);
        WZ=num2str(find(WEIZHI));
        num_digits1 = numel(WZ);
        if num_digits1==num_digits
            cucu=[cucu (10^(num_digits)*str2num(WZ))];
        else
             cucu=[cucu (10^(num_digits-num_digits1)*str2num(WZ))];
        end        
    end
    baocun_mingcheng(i)={cucu};
end

for i=1:length(baocun_mingcheng)
    data=[];
    data=baocun_mingcheng{1,i};
    cc=[];
    data_n=[];data_m=[];
    [data_n,data_m]=sort(data);
    for j=1:length(data)
        cc=[cc num2str(data_n(j))];
    end
    baocun_mingcheng1(i)={cc};
end
    
jilu=[];
for i=1:length(baocun_mingcheng)
    zancun=[];
    zancun=baocun_mingcheng1{1,i};
    chongfu=[];
    chongfu=ismember(baocun_mingcheng1,zancun);
    a=[];
    a=find(chongfu);
    if length(a)>1
        jilu(length(jilu)+1:length(jilu)+length(a)-1)=a(2:end);
    end
end
houxuan(unique(jilu),:)=[];

%Ddd the constant terms
L=size(houxuan,1);
houxuan{L+1,1}='c0';
houxuan{L+1,2}=ones(1,n);
 end