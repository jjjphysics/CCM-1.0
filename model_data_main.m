clc
clear all
lk=linspace(log10(0.001),log10(10),15);
lk=10.^(lk);
count1=0;the_last=[];
for eta=lk
    count1=count1+1;
    n=10000;
    X=zeros(n,2);
    lp=linspace(log10(0.001),log10(3*10^(-2)),40);
    lp=10.^(lp);
    final=[];count=0;
    for epsilong=lp
        count=count+1;
        tempp=[];
        for redo=1:10
            X(1,1)=0.4;X(1,2)=0.2;
            for i=2:n
                X(i,1)=X(i-1,1)*(3.8-3.8*X(i-1,1)-0.05*X(i-1,2));
                X(i,2)=X(i-1,2)*(3.5-3.5*X(i-1,2)-0.1*X(i-1,1));
            end
            si=size(X);
            T=epsilong*randn(si(1,1),2);
            X(:,1)=X(:,1)+eta*T(:,1); X(:,2)=X(:,2)+T(:,2);
            [rou,lpk]=ccm(X(1:si(1,1),1),X(1:si(1,1),2),2,1,1000,2);
            tempp(:,redo,1)=rou(:,1);tempp(:,redo,2)=rou(:,2);
        end
        final(count,1)=epsilong;
        final(count,2)=mean(reshape(tempp(:,:,1),size(tempp,1)*size(tempp,2),1));
        final(count,3)=mean(reshape(tempp(:,:,2),size(tempp,1)*size(tempp,2),1));
        final(count,4)=std(reshape(tempp(:,:,1),size(tempp,1)*size(tempp,2),1));
        final(count,5)=std(reshape(tempp(:,:,2),size(tempp,1)*size(tempp,2),1));
    end
    the_last(:,count1)=final(:,3)-final(:,2);
end