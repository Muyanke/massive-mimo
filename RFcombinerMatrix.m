function Wk=RFcombinerMatrix(N_MS,M_MS,Hk)
% Wk=RFcombinerMatrix(N_MS,M_MS,N_BS,Nc,Np,Xbs,Ybs,Xms,Yms)

W=zeros(N_MS,N_MS);%存储所有的d(w)
W_norm=zeros(N_MS,1);%存储(norm(d'*Hk,1)).^2计算出来的范数
Wk=zeros(N_MS,M_MS);%初始化合成器矩阵Wk，为N_MS*M_MS维度
%Hk=ChannelMatrix(Nc,Np,N_BS,N_MS,Xbs,Ybs,Xms,Yms);
d=zeros(N_MS,1);
for w=1:1:N_MS    
    d=d_w(N_MS,2*pi*(w-1)/N_MS);
    W_norm(w)=(norm(d'*Hk,1)).^2;
    W(:,w)=d;
end
[W_sorted,index]=sort(W_norm,'descend');%降序排列
chosen=index(1:1:M_MS,:);
Wk=1/sqrt(N_MS)*W(:,chosen);%选择前M_MS个Wk