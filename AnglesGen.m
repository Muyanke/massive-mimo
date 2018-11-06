function Angles=AnglesGen(max,min,sigma)
%% 方向角函数，利用截断拉普拉斯分布（分布函数）生成AOA/AOD
%sigma laplacian分布标准差
%max 均值均匀分布范围的最大值
%min 均值均匀分布范围的最小值
%% begin
miu=unifrnd(min,max);%miu 均值
b=sigma/sqrt(2);    %b 尺度参数
a=rand()-0.5;    %生成(-0.5,0.5)区间内均匀分布的随机数;

Angles=miu-b*sign(a).*log(1-2*abs(a)); %生成符合拉普拉斯分布的随机数

%Angles =1/(2*b)*exp(-abs(x-miu)/b);  %此处生成的是概率密度函数
% if(x<=miu)
%     Angles=0.5*exp(-abs(miu-x)/b);
% else
%     Angles=1-0.5*exp(-abs(x-miu)/b);
% end                                   %此处生成的是概率分布函数
%plot(x, Angles)可视化