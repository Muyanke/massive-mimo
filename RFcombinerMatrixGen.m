%% ����RF��ϳɾ���
N_BS = 8;%��վ��������
N_MS = 8;%���ն�������
Ns = 2;%ÿ�����ն˴�����������
%MS��

Nc = 8;%Number of clusters��Ⱥ��
Np = 10;%Number of rays in each clusterÿ����Ⱥ·����


%��������Ԫ����,Y��������Ԫ��Ϊ1����ΪULA����
%�����ϵʽXbs*Ybs=N_BS
%         Xms*Yms=N_MS
Xbs=8;%���Ͷ�����x����Ԫ��
Ybs=1;%���Ͷ�����y����Ԫ��
Xms=8;%���ն�����x����Ԫ��
Yms=1;%���ն�����y����Ԫ��
Dbs=1/2;%��Ԫ���
Dms=1/2;

%����Ǻ����Ǽ������
sigma = 7.5/180*pi;    %***ASSUME the spread are EQUAL.***
sigmaPhi_bs = sigma;   %in azimuth at BS
sigmaPhi_ms = sigma;   %in azimuth at MS
sigmaTheta_bs = sigma; %in elevation at BS
sigmaTheta_ms = sigma; %in elevation at MS

maxPhi_bs = pi/3;  %Azimuth range -60 to 60
maxPhi_ms = pi/3;
maxTheta_bs = pi/9;%Elevation range -20 to 20
maxTheta_ms = pi/9;

%% begin
M_MS=2; % Number of RF chains in Mobilestation
% Wk=zeros(N_MS,N_MS);
% W=zeros(N_MS,1);
% 
Hk=ChannelMatrix(Nc,Np,N_BS,N_MS,Xbs,Ybs,Xms,Yms);
% d=zeros(N_MS,1);
% for w=1:1:N_MS    
%     d=d_w(N_MS,2*pi*(w-1)/N_MS);
%     W(w)=(norm(d'*Hk,1)).^2;
%     Wk(:,w)=d;
% end
% [W_sorted,index]=sort(W,'descend');%��������
% chosen=index(1:1:M_MS,:);
% Wk=1/sqrt(N_MS)*Wk(:,chosen);%ѡ��ǰM_MS��Wk


W=zeros(N_MS,N_MS);
W_norm=zeros(N_MS,1);
Wk=zeros(N_MS,M_MS);
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