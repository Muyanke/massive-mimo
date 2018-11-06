function Angles=AnglesGen(max,min,sigma)
%% ����Ǻ��������ýض�������˹�ֲ����ֲ�����������AOA/AOD
%sigma laplacian�ֲ���׼��
%max ��ֵ���ȷֲ���Χ�����ֵ
%min ��ֵ���ȷֲ���Χ����Сֵ
%% begin
miu=unifrnd(min,max);%miu ��ֵ
b=sigma/sqrt(2);    %b �߶Ȳ���
a=rand()-0.5;    %����(-0.5,0.5)�����ھ��ȷֲ��������;

Angles=miu-b*sign(a).*log(1-2*abs(a)); %���ɷ���������˹�ֲ��������

%Angles =1/(2*b)*exp(-abs(x-miu)/b);  %�˴����ɵ��Ǹ����ܶȺ���
% if(x<=miu)
%     Angles=0.5*exp(-abs(miu-x)/b);
% else
%     Angles=1-0.5*exp(-abs(x-miu)/b);
% end                                   %�˴����ɵ��Ǹ��ʷֲ�����
%plot(x, Angles)���ӻ�