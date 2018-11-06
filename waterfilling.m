clear
clc 
nt=[2,4,8];          %发射天线的数目  
nr=[2,4,8];          %接收天线的数目  
Pt_db=0:5:30;        %信噪比（单位dB)  
PT=10.^(Pt_db/10);   %信噪比的单位转换  
for i=1:3  
    for j=1:length(Pt_db)  
         Cn=[];  
        for k=1:1000                      %仿真时的抽样数量  
            H=(randn(nr(i),nt(i))+sqrt(-1)*randn(nr(i),nt(i)))/sqrt(2);  %瑞利衰落信道矩阵   
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
    ylabel('容量，单位为bps');  
    title('MIMO信道容量及注水算法');  
end  