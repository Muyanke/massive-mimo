clear
clc

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

N_BS = 8;%��վ��������
N_MS = 8;%���ն�������
M_BS=8; %��վRFchain    %M_BS=K*M_MS
M_MS=4; %���ն�RFchain
Ns = 2;%ÿ�����ն˴�����������
K_all = 4; %MS��


User=1:K_all;%��ʼ���û�����Ϊȫ�û�
H_all_k=zeros(N_MS,N_BS,K_all);

%���������û����ŵ�����
for k=1:K_all
    rng(k);%�������������
    H_all_k(:,:,k)=ChannelMatrix(Nc,Np,N_BS,N_MS,Xbs,Ybs,Xms,Yms);
end



while(size(User,2)>1)


K=size(User,2);
M_BS=M_MS*K;     %�û����ϸı���û���RFchain��Ҳ����ı䣬��Ȼ�ᳬ��ά������
H=zeros(N_MS,N_BS,K);

for k=1:K
    H(:,:,k)=H_all_k(:,:,User(k));  %�ı����û��������û���Ӧ���ŵ�����
end

%����PrecoderAndCombinerMatrix��������Ԥ�������ͺϳɾ���֮���Sigma����עˮ
[W,F,B,M,Sigma]=PrecoderAndCombinerMatrix(Nc,Np,N_BS,N_MS,M_MS,M_BS,Ns,Xbs,Ybs,Xms,Yms,K,H);



Sigma_n=[];
for k=1:K
    Sig=diag(Sigma(:,:,k));
    Sigma_n=[Sigma_n;Sig];
end

    
%% עˮ�㷨���ʷ��� 
% û�п��Ǽ�Ȩϵ��w_k   

SNR=0;%��λdB  SNR=P/sigma^2;
Sigma_n=(10^(SNR./10))./(K*Ns)*Sigma_n;

%�û�����Ȩ��
% rng(K+1)
% weight_k=rand(K,1);

%�û�֮���Ȩ��
%����cvx��עˮ�㷨�õ����ʷ������
weight_k=ones(K,1);
tmp=repmat(weight_k',Ns,1);
weight_n=reshape(tmp,K*Ns,1);


cvx_begin
    variable Lambda_n(K*Ns) 
    minimize(-sum(weight_n.*log(1+Lambda_n.*Sigma_n)/log(2)))
    subject to
        sum(Lambda_n)==K*Ns;
        Lambda_n>=0;
cvx_end
cvx_obj = cvx_optval;

Gama=Lambda_n;


%�õ����еĹ��ʷ������
Gama_k=zeros(Ns,Ns,K);
Rate_k=zeros(K,1);

for k=1:K
    Gama_k(:,:, k)=diag(Gama(2*k-1:2*k,1));
    Rate_k(k,1)=log2(det(eye(Ns)+((10^(SNR./10))./(K*Ns)).*Gama_k(:,:, k)*Sigma(:,:,k).^2));
end



%% �ڵ㷨 

% P_k_max=0.5;    %��λW
% Power_n=((10^(SNR./10))./(K*Ns))*Gama;%ÿ���������Ϸ���Ĺ���
% R_k_min=5;     %��λbps/Hz
% TempRk=repmat(Rate_k',Ns,1);
% Rate_k_n=reshape(TempRk,K*Ns,1);
% cvx_begin
%     variable alpha_k(K*Ns) 
%     minimize(sum((alpha_k-1).^2))
%     subject to
%         Power_n<=P_k_max;
%         alpha_k.^2.*R_k_min-Rate_k_n<=0;
% cvx_end
% cvx_obj = cvx_optval;



P_k_max=0.5;    %��λW

Power_n=((10^(SNR./10))./(K*Ns))*Gama;%ÿ���������Ϸ���Ĺ���
Power_k=zeros(K,1);
for i=1:K
    Power_k(i)=sum(Power_n(i:i*Ns)); %ÿ���û�����Ĺ���
end

R_k_min=5;     %��λbps/Hz


cvx_begin
    variable alpha_k(K) 
    minimize(sum((alpha_k-1).^2))
    subject to
        Power_k<=P_k_max;
        alpha_k.^2.*R_k_min-Rate_k<=0;
cvx_end
cvx_obj = cvx_optval;

if alpha_k==1
    break;
else
    [min_alpha,index]=min(alpha_k);
    User(index)=[];
end
end
User
Power_k