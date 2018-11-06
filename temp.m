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
M_BS=8; %基站RFchain    %M_BS=K*M_MS
M_MS=4; %接收端RFchain
Ns = 2;%每个接收端处理数据流数
K_all = 4; %MS数


User=1:K_all;%初始化用户集合为全用户
H_all_k=zeros(N_MS,N_BS,K_all);

%生成所有用户的信道矩阵
for k=1:K_all
    rng(k);%设置随机数种子
    H_all_k(:,:,k)=ChannelMatrix(Nc,Np,N_BS,N_MS,Xbs,Ybs,Xms,Yms);
end



while(size(User,2)>1)


K=size(User,2);
M_BS=M_MS*K;     %用户集合改变后用户端RFchain数也必须改变，不然会超过维度索引
H=zeros(N_MS,N_BS,K);

for k=1:K
    H(:,:,k)=H_all_k(:,:,User(k));  %改变后的用户集合中用户对应的信道矩阵
end

%调用PrecoderAndCombinerMatrix函数计算预编码矩阵和合成矩阵，之后对Sigma进行注水
[W,F,B,M,Sigma]=PrecoderAndCombinerMatrix(Nc,Np,N_BS,N_MS,M_MS,M_BS,Ns,Xbs,Ybs,Xms,Yms,K,H);



Sigma_n=[];
for k=1:K
    Sig=diag(Sigma(:,:,k));
    Sigma_n=[Sigma_n;Sig];
end

    
%% 注水算法功率分配 
% 没有考虑加权系数w_k   

SNR=0;%单位dB  SNR=P/sigma^2;
Sigma_n=(10^(SNR./10))./(K*Ns)*Sigma_n;

%用户不等权重
% rng(K+1)
% weight_k=rand(K,1);

%用户之间等权重
%利用cvx解注水算法得到功率发配矩阵
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


%得到所有的功率分配矩阵
Gama_k=zeros(Ns,Ns,K);
Rate_k=zeros(K,1);

for k=1:K
    Gama_k(:,:, k)=diag(Gama(2*k-1:2*k,1));
    Rate_k(k,1)=log2(det(eye(Ns)+((10^(SNR./10))./(K*Ns)).*Gama_k(:,:, k)*Sigma(:,:,k).^2));
end



%% 内点法 

% P_k_max=0.5;    %单位W
% Power_n=((10^(SNR./10))./(K*Ns))*Gama;%每个数据流上分配的功率
% R_k_min=5;     %单位bps/Hz
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



P_k_max=0.5;    %单位W

Power_n=((10^(SNR./10))./(K*Ns))*Gama;%每个数据流上分配的功率
Power_k=zeros(K,1);
for i=1:K
    Power_k(i)=sum(Power_n(i:i*Ns)); %每个用户分配的功率
end

R_k_min=5;     %单位bps/Hz


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