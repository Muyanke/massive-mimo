function a_ULA=ArrayResponseVec(N,M,D,Phi,Theta)
%% %function ArrayResponseVec
%UPA阵列响应向量
%N X轴天线数
%M Y轴天线数，M=1时退化为ULA
%D 阵元间距
%Phi 方向角
%Theta 仰角
%% begin
a=zeros(N*M,1);
for n=1:N
    for m=1:M
    a(n*m,1)=exp(1j*2*pi*D*((m-1)*sin(Theta)*cos(Phi)+ (n-1)*sin(Theta)*sin(Phi)));
    end
end
a_ULA=(1/sqrt(N*M))*a;
