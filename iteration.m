clear
clc

Nc = 8;%Number of clusters��Ⱥ��
Np = 10;%Number of rays in each clusterÿ����Ⱥ·����


%��������Ԫ����,Y��������Ԫ��Ϊ1����ΪULA����
%�����ϵʽ Xbs*Ybs=N_BS
%          Xms*Yms=N_MS
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
M_BS=8; %��վRFchain
M_MS=2; %���ն�RFchain
Ns = 2;%ÿ�����ն˴�����������

K = 4; %MS��%M_BS=K*M_MS






%% RF����
%����RF��Ԥ�������ͺϳ�������
H=zeros(N_MS,N_BS,K);
W=zeros(N_MS,M_MS,K);
H_int=[];
for k=1:K
    rng(k);%������������ӣ�������
    H(:,:,k)=ChannelMatrix(Nc,Np,N_BS,N_MS,Xbs,Ybs,Xms,Yms);
    W(:,:,k)=RFcombinerMatrix(N_MS,M_MS,H(:,:,k));
    H_int=[H_int;W(:,:,k)'*H(:,:,k)];%Ӧ��Ϊ����ƴ��K��  
end
angH_int=angle(H_int');
F=1/sqrt(N_BS)*exp(1j*angH_int);

%Heq=H_int*F;%��Ч�����ŵ�


 

%% ��������

%����������Ԥ�������B�ͺϳ�������M

    B=zeros(M_BS,Ns,K);%M_BS=K*M_MS������Heq��һ������
    M=zeros(M_MS,Ns,K);
    
    H_withoutk=cell(1,K);
    Beq=cell(1,K);
    Sigma=zeros(Ns,Ns,K);
    for k=1:K
        for j=1:K
            H_withoutk{1,j}=(W(:,:,j)'*H(:,:,j)*F)';
            if j==k
            H_withoutk{1,j}=[];
            end
        end
    H_without_k=(cell2mat(H_withoutk))';%����Hk~,��ȥk��(K-1)*M_MS�У�M_BS��
   
    %������SVD�ֽ�
    
    r=rank(H_without_k);
    [U,S,V]=svd(H_without_k);
    V1=V(:,1:r);
    V2=V(:,r+1:M_BS); %M_BS=K*M_MS
    
    %��һ��svd
    Hk_wave=W(:,:,k)'*H(:,:,k)*F;
    [P,E,Q]=svd(Hk_wave*V2);
    P1=P(:,1:Ns);
    Q1=Q(:,1:Ns);
    M(:,:,k)=P1;
    B(:,:,k)=V2*Q1;
    Beq{1,k}=B(:,:,k);
   %% עˮ�����ɹ��ʷ������
   Sigma(:,:,k)=E(1:Ns,1:Ns);
   
      
    end
B=cell2mat(Beq);


Sigma_n=[];
for i=1:K
    Sig=diag(Sigma(:,:,k));
    Sigma_n=[Sigma_n;Sig];
end

    
%% עˮ�㷨���ʷ��� 
% û�п��Ǽ�Ȩϵ��w_k   

SNR=0;%��λdB  
Sigma_n=(10^(SNR./10))./(K*Ns)*Sigma_n;


% rng(K+1)
% weight_k=rand(K,1);
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
    Power_k(i)=sum(Power_n(i:i*Ns)); 
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






