function [W,F,B,M,Sigma]=PrecoderAndCombinerMatrix(Nc,Np,N_BS,N_MS,M_MS,M_BS,Ns,Xbs,Ybs,Xms,Yms,K,H)
%% 参数说明
% Nc   Number of clusters集群数
% Np   Number of rays in each cluster每个集群路径数


% Xbs   发送端天线x轴阵元数
% Ybs   发送端天线y轴阵元数
% Xms   接收端天线x轴阵元数
% Yms   接收端天线y轴阵元数
% 阵列中阵元个数,Y轴天线阵元数为1，即为ULA阵列
% 满足关系式Xbs*Ybs=N_BS
%          Xms*Yms=N_MS


% N_BS  基站端天线数
% N_MS  接收端天线数
% M_BS  基站RFchain
% M_MS  接收端RFchain

% K     MS数
% Ns    每个接收端处理数据流数

%% RF域处理
%生成RF域预编码矩阵和合成器矩阵

W=zeros(N_MS,M_MS,K);
H_int=[];
for k=1:K
    
    W(:,:,k)=RFcombinerMatrix(N_MS,M_MS,H(:,:,k));
    H_int=[H_int;W(:,:,k)'*H(:,:,k)];%应该为纵向拼接K行  
end
angH_int=angle(H_int');
F=1/sqrt(N_BS)*exp(1j*angH_int);

%Heq=H_int*F;%等效基带信道


 

%% 数字域处理

%生成数字域预编码矩阵B和合成器矩阵M

    B=zeros(M_BS,Ns,K);%M_BS=K*M_MS，导致Heq是一个方阵
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
    H_without_k=(cell2mat(H_withoutk))';  %生成Hk~,除去k，(K-1)*M_MS行，M_BS列
   
    %下面做SVD分解
    
    r=rank(H_without_k);
    [U,S,V]=svd(H_without_k);
    %V1=V(:,1:r);
    V2=V(:,r+1:M_BS); %M_BS=K*M_MS
    
    %进一步svd
    Hk_wave=W(:,:,k)'*H(:,:,k)*F;
    [P,E,Q]=svd(Hk_wave*V2);
    P1=P(:,1:Ns);
    Q1=Q(:,1:Ns);
    M(:,:,k)=P1;
    B(:,:,k)=V2*Q1;
    %Beq{1,k}=B(:,:,k);
   
    %% 注水法生成功率分配矩阵需要的Sigma矩阵
   Sigma(:,:,k)=E(1:Ns,1:Ns);
   
      
    end
   %B=cell2mat(Beq);