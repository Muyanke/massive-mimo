%数字域预编码矩阵和合成矩阵生成，用
clear
clc

Nc = 8;%Number of clusters集群数
Np = 10;%Number of rays in each cluster每个集群路径数


%阵列中阵元个数,Y轴天线阵元数为1，即为ULA阵列
%满足关系式Xbs*Ybs=N_BS
%         Xms*Yms=N_MS
Xbs=8;%发送端天线x轴阵元数
Ybs=1;%发送端天线y轴阵元数
Xms=8;%接收端天线x轴阵元数
Yms=1;%接收端天线y轴阵元数
Dbs=1/2;%阵元间距
Dms=1/2;

%方向角和仰角假设相等
sigma = 7.5/180*pi;    %***ASSUME the spread are EQUAL.***
sigmaPhi_bs = sigma;   %in azimuth at BS
sigmaPhi_ms = sigma;   %in azimuth at MS
sigmaTheta_bs = sigma; %in elevation at BS
sigmaTheta_ms = sigma; %in elevation at MS

maxPhi_bs = pi/3;  %Azimuth range -60 to 60
maxPhi_ms = pi/3;
maxTheta_bs = pi/9;%Elevation range -20 to 20
maxTheta_ms = pi/9;

N_BS = 8;%基站端天线数
N_MS = 8;%接收端天线数
M_BS=8; %基站RFchain
M_MS=2; %接收端RFchain
Ns = 2;%每个接收端处理数据流数
%MS数
K = 4; 




%H_int=Wk'*Hk;
% angH_int=angle(H_int);
% Fk=1/sqrt(N_BS)*exp(angH_int);

%F=zeros(N_BS,K*M_MS);
%% RF域处理
%生成RF域预编码矩阵和合成器矩阵
H=zeros(N_MS,N_BS,K);
W=zeros(N_MS,M_MS,K);
H_int=[];
for k=1:K
    H(:,:,k)=ChannelMatrix(Nc,Np,N_BS,N_MS,Xbs,Ybs,Xms,Yms);
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
    for k=1:K
        
        for j=1:K
            H_withoutk{1,j}=(W(:,:,j)'*H(:,:,j)*F)';
            if j==k
            H_withoutk{1,j}=[];
            end
        end
    H_without_k=(cell2mat(H_withoutk))';%生成Hk~,除去k，(K-1)*M_MS行，M_BS列
   
    %下面做SVD分解
    
    r=rank(H_without_k);
    [U,S,V]=svd(H_without_k);
    V1=V(:,1:r);
    V2=V(:,r+1:M_BS); %M_BS=K*M_MS
    
    %进一步svd
    Hk_wave=W(:,:,k)'*H(:,:,k)*F;
    [P,E,Q]=svd(Hk_wave*V2);
    P1=P(:,1:Ns);
    Q1=Q(:,1:Ns);
    
    M(:,:,k)=P1;
    B(:,:,k)=V2*Q1;
    Beq{1,k}=B(:,:,k);
    
    %% 注水法生成功率分配矩阵，waterfilling power allocation
    %Gama=B;
    
    
    end
    B=cell2mat(Beq);