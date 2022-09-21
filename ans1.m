close all
clc;
clear;

global m1 m2 mf cx k 
 f rho g r
m1=4866; % 浮子
mf=1335.535; % 附加垂荡质量
cx=656.3616; % 兴波阻尼系数;
m2=2433; % 振子
k=80000; % 弹簧刚度
omega=1.4005; % 激励频率
f=6250; % 激励力幅值
rho=1025; % 海水密度
g=9.8;
r=1; % 半径
T=2*pi/omega; % 一个波浪周期
Ttotal=40*T; 
dt=0.2; 
tspan=0:dt:Ttotal; 
y0=[0;0;0;0];
[tt1,yy1]=ode45(@myode1,tspan,y0);
[tt2,yy2]=ode45(@myode2,tspan,y0); % 分别求解两种阻尼系数情况下的解

disp('第一种情况，定常阻尼系数:')
for time=[10 20 40 60 100]
    n=time/dt+1;
    d1_fu=yy1(n,1);
    v1_fu=yy1(n,2);
    d1_zh=yy1(n,3);
    v1_zh=yy1(n,4); % 第一种情况下浮子和振子位移和速度在各个时间点处的值
    
    fprintf('\t时间为%ds时，浮子位移为%f/m,速度为%f/(m/s);振子位移为%f/m,速度为%f/(m/s).\n',time,d1_fu,v1_fu,d1_zh,v1_zh)
end

disp('第二种情况，非定常阻尼系数:')
for time=[10 20 40 60 100]
    n=time/dt+1;
    d1_fu=yy2(n,1);
    v1_fu=yy2(n,2);
    d1_zh=yy2(n,3);
    v1_zh=yy2(n,4); % 第二种情况下浮子和振子位移和速度在各个时间点处的值
    
    fprintf('\t时间为%ds时，浮子位移为%f/m,速度为%f/(m/s);振子位移为%f/m,速度为%f/(m/s).\n',time,d1_fu,v1_fu,d1_zh,v1_zh)
end
%%
function dy=myode1(t,y)
global m1 m2 mf cx k omega f rho g r
c=10000; % （1） 常数阻尼系数
dy=zeros(4,1);
dy(1)=y(2);
dy(2)=-cx/(m1+mf)*y(2)-pi*r^2*rho*g/(m1+mf)*y(1)-c/(m1+mf)*(y(2)-y(4))-k/(m1+mf)*(y(1)-y(3))+f/(m1+mf)*cos(omega*t);
dy(3)=y(4);
dy(4)=-c/m2*(y(4)-y(2))-k/m2*(y(3)-y(1));
end

function dy=myode2(t,y)
global m1 m2 mf cx k omega f rho g r
c=10000*sqrt(abs(y(2)-y(4))); % （2） 非定常阻尼系数
dy=zeros(4,1);
dy(1)=y(2);
dy(2)=-cx/(m1+mf)*y(2)-pi*r^2*rho*g/(m1+mf)*y(1)-c/(m1+mf)*(y(2)-y(4))-k/(m1+mf)*(y(1)-y(3))+f/(m1+mf)*cos(omega*t);
dy(3)=y(4);
dy(4)=-c/m2*(y(4)-y(2))-k/m2*(y(3)-y(1));
end
