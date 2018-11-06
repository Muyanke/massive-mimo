%% 生成RF域合成矩阵
N_BS = 8;%基站端天线数
N_MS = 8;%接收端天线数
Ns = 2;%每个接收端处理数据流数
%MS数

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
% [W_sorted,index]=sort(W,'descend');%降序排列
% chosen=index(1:1:M_MS,:);
% Wk=1/sqrt(N_MS)*Wk(:,chosen);%选择前M_MS个Wk


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
[W_sorted,index]=sort(W_norm,'descend');%降序排列
chosen=index(1:1:M_MS,:);
Wk=1/sqrt(N_MS)*W(:,chosen);%选择前M_MS个Wk