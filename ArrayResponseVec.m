function a_ULA=ArrayResponseVec(N,M,D,Phi,Theta)
%% %function ArrayResponseVec
%UPA������Ӧ����
%N X��������
%M Y����������M=1ʱ�˻�ΪULA
%D ��Ԫ���
%Phi �����
%Theta ����
%% begin
a=zeros(N*M,1);
for n=1:N
    for m=1:M
    a(n*m,1)=exp(1j*2*pi*D*((m-1)*sin(Theta)*cos(Phi)+ (n-1)*sin(Theta)*sin(Phi)));
    end
end
a_ULA=(1/sqrt(N*M))*a;
