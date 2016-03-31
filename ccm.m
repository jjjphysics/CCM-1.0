function [final_rou,lpk,record,re_1]=ccm(X1,X2,X3,X4,X5,X6)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ccm calculates the directed dynamical infulence between two time series.
%    
%     [final_rou,lpk,recod,re_1]=ccm(X1,X2,X3,X4,X5,X6)
%
% Inputs:
%         X1:  data sets for first sample(n*1)
%         X2:  data sets for seconds sample(n*1)
%         X3:  the embeding dimensions(E)
%         X4:  the time delay interval(tau)
%         X5:  the largest time series used(L) 
%         X6:  the how much point you want to calculate between L(A integer less than L)
%         the time sets of X1 and X2 must have the same length
%  Outputs:
%         final_rou: the final result of correlation
%         final_rou(i,1): the value \rou_{Y|M_X}
%         final_rou(i,2): the value \rou_{X|M_Y}
%         lpk: the length of time series used calculate ccm
%         record: time-delay lagged -coordinate vector series
%         re_1: Euclidean distance between any two point in lagged-coordinate embedding manifold

%
%  CCM1.0
%  by Junjie Jiang, Arizona State University, Tempe, AZ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sit1=size(X1);
if X5>sit1(1,1)
    disp('L is longer than data sets!L must short than data sets.');
end
X(:,1)=X1;X(:,2)=X2;
final_rou=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol=X4;
count=0;
% lpk=linspace(log10(3*X3),log10(X5),X6);
% lpk=10.^(lpk);
% lpk=floor(lpk);
lpk=1000;
for j=lpk
    count=count+1;
    for re_do=1:10
    %% find random initial point and get two time-delay lagged-coordinate vector series
        
        q=floor(X3*tol+(X5-j)*rand);
        record=[];
        for temp1=1:X3
            select1=q-tol*(temp1-1);
            record(:,temp1)=X(select1:j+select1,1);
            record(:,temp1+X3)=X(select1:j+select1,2);
        end
   %% calculate Euclidean distance between any two point in lagged-coordinate embedding manifold 
    re_1=[];    
    for kl=1:2
        for temp1=1:j+1
            for temp2=1:j+1
                 value=0;
                 for temp3=1:X3
                    value=value+(record(temp1,temp3+X3*(kl-1))-record(temp2,temp3+X3*(kl-1)))^2;   
                 end
                 re_1(temp2,temp1,kl)=sqrt(value);             
            end
        end
    end
    
    %% combine lagged-coordinate vector series with Euclidean distance matrix  
    sig=size(re_1);
    for jk=1:2
        for lk=1:2*X3
            re_1(:,sig(1,2)+lk,jk)=record(:,lk);
        end
    end
    %% Cross mapping calculate(E+1 nearest neighbors finding, lagged-coordinate vectors estimate and correlation between estimate lagged-coordinate vectors and lagged-coordinate vectors)     
    for temp1=1:2
        number=[];
        wq=j-X3-1;
        for w=(X3+2):j
            sit=size(re_1);
            temp3=re_1(:,w,temp1);temp3(:,2:(X3+1))=record(:,(2*X3-X3*temp1+1):(2*X3-X3*(temp1-1)));
            temp3=sortrows(temp3(:,:),1);
            temp4=temp3(2:X3+2,:);
            sum1=0;temp_1=[];
            for po=1:X3+1
                sum1=sum1+exp(-temp4(po,1)/temp4(1,1));
                temp_1(po,1)=exp(-temp4(po,1)/temp4(1,1));
            end
            temp_1=temp_1/sum1;
            temp_1=temp_1';
            temp_2=temp4(1:X3+1,2:(X3+1));
            y_bar=temp_1*temp_2;
            number(w-X3-1,:)=y_bar;
            result(w-X3-1,1,temp1)=y_bar(1,1);result(w-X3-1,2,temp1)=y_bar(1,2);
        end
        number(:,(X3+1):2*X3)=re_1((X3+2):j,(sit(1,2)-X3*(temp1-1))-X3+1:(sit(1,2)-X3*(temp1-1)));
        sum2=zeros(X3,5);
        for qr=1:wq
            for qt=1:X3
                sum2(qt,1)=sum2(qt,1)+number(qr,qt);
                sum2(qt,2)=sum2(qt,2)+number(qr,qt+X3);
                sum2(qt,3)=sum2(qt,3)+(number(qr,qt))^2;
                sum2(qt,4)=sum2(qt,4)+(number(qr,qt+X3))^2;
                sum2(qt,5)=sum2(qt,5)+number(qr,qt)*number(qr,qt+X3);
            end
        end
        rou1=0;
        for qy=1:X3
            rou1=rou1+((sum2(qy,5)/(wq))-(sum2(qy,1)/(wq))*(sum2(qy,2)/(wq)))/((sqrt((sum2(qy,3)/(wq))-(sum2(qy,1)/(wq))^2))*(sqrt((sum2(qy,4)/(wq))-(sum2(qy,2)/(wq))^2)));
        end
        rou(1,temp1)=rou1/X3;                    
    end
%     sum_temp1(re_do,1)=rou(1,1);
%     sum_temp2(re_do,1)=rou(1,2);
    final_rou(re_do,1)=rou(1,1);
    final_rou(re_do,2)=rou(1,2);
    end
%     final(count,1)=mean(sum_temp1,1);
%     final(count,2)=std(sum_temp1,1,1);
%     final(count,3)=mean(sum_temp2,1);
%     final(count,4)=std(sum_temp2,1,1);
end
end
