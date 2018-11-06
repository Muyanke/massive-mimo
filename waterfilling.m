clear
clc 
nt=[2,4,8];          %�������ߵ���Ŀ  
nr=[2,4,8];          %�������ߵ���Ŀ  
Pt_db=0:5:30;        %����ȣ���λdB)  
PT=10.^(Pt_db/10);   %����ȵĵ�λת��  
for i=1:3  
    for j=1:length(Pt_db)  
         Cn=[];  
        for k=1:1000                      %����ʱ�ĳ�������  
            H=(randn(nr(i),nt(i))+sqrt(-1)*randn(nr(i),nt(i)))/sqrt(2);  %����˥���ŵ�����   
            [U,D,V]=svd(H);  
            A=D'*D;  
            r=rank(H);  
            a=1./diag(A);  
            s=0;  
            b=sort(a);  
            for m=1:(r-1)  
                s=s+m*(b(m+1)-b(m));  
               if s>PT(j)  
                   v=b(m+1)-(s-PT(j))/m;  
                   break;  
                end  
            end  
            if s<PT(j)  
               v=b(r)+(PT(j)-s)/r;  
            end   
            for n=1:r  
               x(n)=max(v-a(n),0);  
            end  
            x=x(1:r);  
            X=diag(x);  
            Rx=V*X*V';  
            I=diag(ones(1,nr(i)));  
            c=log(det(I+H*Rx*H'));  
            Cn=[Cn,c];             
        end  
        y(j)=real(sum(Cn)/1000);  
    end  
    plot(Pt_db,y);  
    hold on;  
    xlabel('Pt_db');  
    ylabel('��������λΪbps');  
    title('MIMO�ŵ�������עˮ�㷨');  
end  