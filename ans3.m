close all
clc;
clear;
global m1 m2 mf cx k omega f rho g r J1 J_fu C_hui C_xing k_niu k_zu L
m1=4866; % 浮子
mf=1028.876; % 附加垂荡质量 0
J1=compute_J; % 计算得到的浮子绕质心的转动惯量!0
J_fu=7001.914; % 附加转动惯量;0
cx=683.4558; % 兴波阻尼系数;0
C_xing=654.3383; % 摇荡兴波阻尼系数 0
C_hui=8890.7; % 静水恢复力矩系数 0
m2=2433; % 振子
k=80000; % 弹簧刚度 0
k_niu=250000; % 扭转弹簧刚度
k_zu=1000; % 旋转阻尼器的阻尼系数 0
omega=1.7152; % 激励频率 0
f=3640; % 激励力幅值
L=1690; % 摇荡激励力矩幅值 0
rho=1025; % 海水密度
g=9.8;
r=1; % 半径
T=2*pi/omega; % 一个波浪周期 0
Ttotal=40*T; % 40个周期
dt=0.2; % 间隔0.2s
tspan=0:dt:Ttotal; % 总的求解时间区间 0
y0=[0;0;0;0;0;0;0;0]; % 初始静止！
[tt1,yy1]=ode45(@myode1,tspan,y0); % 求解! 0
%%% 其中，     tt:各个时间点 0
%%% yy(:,1)第一列: 浮子位移
%%% yy(:,2)第二列：浮子速度 0
%%% yy(:,3)第三列: 振子位移
%%% yy(:,4)第四列：振子速度
%%% yy(:,5)第5列: 浮子角位移 0
%%% yy(:,6)第6列：浮子角速度
%%% yy(:,7)第7列: 振子角位移 0
%%% yy(:,8)第8列：振子角速度
figure
set(gcf,'Position',[100 100 700 200]) 
plot(tt1,yy1(:,1),tt1,yy1(:,3)); % 绘制浮子和振子位移 0
xlabel('time/s')
ylabel('displacement/m')
legend('浮子','振子')
title('浮子与振子垂荡位移')
figure
set(gcf,'Position',[100 400 700 400])
subplot(211)
plot(tt1,yy1(:,5),tt1,yy1(:,7)); % 绘制浮子和振子角位移 0
xlabel('      time/s')
ylabel('angle/rad')
legend('浮子','振子')
title('浮子与振子摇荡角位移')
subplot(212)
plot(tt1,yy1(:,5)-yy1(:,7),'r'); % 绘制浮子和振子角位移之差 0
xlabel('time/s')
ylabel('angle/rad')
title('浮子与振子摇荡角位移之差值')

disp('浮子与振子垂荡位移')
for time=[10 20 40 60 100]
    n=time/dt+1;
    d_fu=yy1(n,1);
    v_fu=yy1(n,2);
    d_zh=yy1(n,3);
    v_zh=yy1(n,4); % 浮子和振子位移和速度在各个时间点处的值 0
    angle_fu=yy1(n,5);
    anglev_fu=yy1(n,6);
    angle_zh=yy1(n,7);
    anglev_zh=yy1(n,8); % 浮子和振子角位移和角速度在各个时间点处的值 0
    fprintf('\t时间为%ds时，浮子位移为%f/m,速度为%f/(m/s);振子位移为%f/m,速度为%f/(m/s);\n\t\t浮子角位移为%f/rad,角速度为%f/(rad/s);振子角位移为%f/rad,角速度为%f/(rad/s).\n',.....
        time,d_fu,v_fu,d_zh,v_zh,angle_fu,anglev_fu,angle_zh,anglev_zh);
end



%% 定义odes
function dy=myode1(t,y)
global m1 m2 mf cx k omega f rho g r J1 J_fu C_hui C_xing k_niu k_zu L
h0=2;
c=10000; % （1） 常数阻尼系数 0
dy=zeros(8,1);
dy(1)=y(2);
dy(2)=-cx/(m1+mf)*y(2)-(pi*r^2*rho*g/(m1+mf)*y(1)*(y(1)>=h0-3)+pi*r^2*rho*g/(m1+mf)*(h0-3)*(y(1)<h0-3))-c/(m1+mf)*(y(2)-y(4))-k/(m1+mf)*(y(1)-y(3))+f/(m1+mf)*cos(omega*t);
dy(3)=y(4);
dy(4)=-c/m2*(y(4)-y(2))-k/m2*(y(3)-y(1));
dy(5)=y(6);
dy(6)=(L*cos(omega*t)-(C_hui*y(5)+C_xing*y(6)+k_niu*(y(5)-y(7))+k_zu*(y(6)-y(8))))/(J1+J_fu); % 关于浮子转角方程 0
dy(7)=y(8);
J2=(0.5-m2*g/k+y(3)-y(1))^2*m2; % 此时，振子关于铰接点的转动惯量 0
dy(8)=-(k_niu*(y(7)-y(5))+k_zu*(y(8)-y(6))+2*m2*(y(4)-y(2))*y(8)*(0.5-m2*g/k+y(3)-y(1)))/J2; % 关于振子的转角方程 0
end

