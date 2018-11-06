function [W,F,B,M,Sigma]=PrecoderAndCombinerMatrix(Nc,Np,N_BS,N_MS,M_MS,M_BS,Ns,Xbs,Ybs,Xms,Yms,K,H)
%% ����˵��
% Nc   Number of clusters��Ⱥ��
% Np   Number of rays in each clusterÿ����Ⱥ·����


% Xbs   ���Ͷ�����x����Ԫ��
% Ybs   ���Ͷ�����y����Ԫ��
% Xms   ���ն�����x����Ԫ��
% Yms   ���ն�����y����Ԫ��
% ��������Ԫ����,Y��������Ԫ��Ϊ1����ΪULA����
% �����ϵʽXbs*Ybs=N_BS
%          Xms*Yms=N_MS


% N_BS  ��վ��������
% N_MS  ���ն�������
% M_BS  ��վRFchain
% M_MS  ���ն�RFchain

% K     MS��
% Ns    ÿ�����ն˴�����������

%% RF����
%����RF��Ԥ�������ͺϳ�������

W=zeros(N_MS,M_MS,K);
H_int=[];
for k=1:K
    
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
    H_without_k=(cell2mat(H_withoutk))';  %����Hk~,��ȥk��(K-1)*M_MS�У�M_BS��
   
    %������SVD�ֽ�
    
    r=rank(H_without_k);
    [U,S,V]=svd(H_without_k);
    %V1=V(:,1:r);
    V2=V(:,r+1:M_BS); %M_BS=K*M_MS
    
    %��һ��svd
    Hk_wave=W(:,:,k)'*H(:,:,k)*F;
    [P,E,Q]=svd(Hk_wave*V2);
    P1=P(:,1:Ns);
    Q1=Q(:,1:Ns);
    M(:,:,k)=P1;
    B(:,:,k)=V2*Q1;
    %Beq{1,k}=B(:,:,k);
   
    %% עˮ�����ɹ��ʷ��������Ҫ��Sigma����
   Sigma(:,:,k)=E(1:Ns,1:Ns);
   
      
    end
   %B=cell2mat(Beq);