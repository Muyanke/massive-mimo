function Hk=ChannelMatrix(Nc,Np,N_BS,N_MS,Xbs,Ybs,Xms,Yms)
%% %ChannelMatrix Generate
%Nc  Number of clusters��Ⱥ��
%Np  Number of rays in each clusterÿ����Ⱥ·����
%N_MS Number of antennas in mobilestation
%N_BS Number of antennas in basestaion
% Xbs���Ͷ�����x����Ԫ��
% Ybs���Ͷ�����y����Ԫ��
% Xms���ն�����x����Ԫ��
% Yms���ն�����y����Ԫ��
% ��������Ԫ����,Y��������Ԫ��Ϊ1����ΪULA����

%% begin
Dbs=1/2;%��Ԫ���
Dms=1/2;

sigma = 7.5/180*pi;%***ASSUME the spread are EQUAL.***
sigmaPhi_bs = sigma;%in azimuth at BS
sigmaPhi_ms = sigma;%in azimuth at MS
sigmaTheta_bs = sigma;%in elevation at BS
sigmaTheta_ms = sigma;%in elevation at MS

maxPhi_bs = pi/3;  %Azimuth range -60 to 60
maxPhi_ms = pi/3;
maxTheta_bs = pi/9;%Elevation range -20 to 20
maxTheta_ms = pi/9;


%% begin
a=zeros(Nc,Np);
% Hkil=zeros(Xms*Yms,Xbs*Ybs,Nc*Np);
% Hk=zeros(Xms*Yms,Xbs*Ybs);
Hkil=zeros(N_MS,N_BS,Nc*Np);
Hk=zeros(N_MS,N_BS);
for i=1:Nc
    for l=1:Np
        a(i,l)=sqrt(1/2)*(randn+1i*randn);   %��i��Ⱥ��l·���ĸ����棬���Ӹ���˹�ֲ�
        Phi_bs=AnglesGen(maxPhi_bs,-maxPhi_bs,sigmaPhi_bs);
        Theta_bs=AnglesGen(maxTheta_bs,-maxTheta_bs,sigmaTheta_bs);
        a_BS=ArrayResponseVec(Xbs,Ybs,Dbs,Phi_bs,Theta_bs);
        Phi_ms=AnglesGen(maxPhi_ms,-maxPhi_ms,sigmaPhi_ms);
        Theta_ms=AnglesGen(maxTheta_ms,-maxTheta_ms,sigmaTheta_ms);
        a_MS=ArrayResponseVec(Xms,Yms,Dms,Phi_ms,Theta_ms);
        Hkil(:,:,i*l)=a(i,l).*a_MS*(a_BS)';
        Hk=Hk+Hkil(:,:,i*l);
    end
end
Hk=sqrt(N_BS*N_MS/(Nc*Np)).*Hk;
