function Wk=RFcombinerMatrix(N_MS,M_MS,Hk)
% Wk=RFcombinerMatrix(N_MS,M_MS,N_BS,Nc,Np,Xbs,Ybs,Xms,Yms)

W=zeros(N_MS,N_MS);%�洢���е�d(w)
W_norm=zeros(N_MS,1);%�洢(norm(d'*Hk,1)).^2��������ķ���
Wk=zeros(N_MS,M_MS);%��ʼ���ϳ�������Wk��ΪN_MS*M_MSά��
%Hk=ChannelMatrix(Nc,Np,N_BS,N_MS,Xbs,Ybs,Xms,Yms);
d=zeros(N_MS,1);
for w=1:1:N_MS    
    d=d_w(N_MS,2*pi*(w-1)/N_MS);
    W_norm(w)=(norm(d'*Hk,1)).^2;
    W(:,w)=d;
end
[W_sorted,index]=sort(W_norm,'descend');%��������
chosen=index(1:1:M_MS,:);
Wk=1/sqrt(N_MS)*W(:,chosen);%ѡ��ǰM_MS��Wk