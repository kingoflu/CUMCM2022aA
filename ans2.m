clc;
clear;
global m1 m2 mf cx k omega f rho g r
m1=4866;
mf=1165.992;%
cx=167.8395;%
m2=2433;
k=80000;
omega=2.2143;
f=4890;%
rho=1025;
g=9.8;
r=1;
T=2*pi/omega;
Ttotal=40*T;
dt=0.1;
tspan=0:dt:Ttotal;
y0=[0;0;0;0];
nc=100;
C_constant_damp=linspace(0,100000,nc);
P1=zeros(size(C_constant_damp));
hw1=waitbar(0,'第一种情况开始计算...');
pause(1)
for i=1:nc
    c=C_constant_damp(i);   
    [tt1,yy1]=ode45(@(t,y)myode1(t,y,c),tspan,y0);
    v_relative=yy1(:,2)-yy1(:,4);
    pl(i)=1/2*c*trapz(tt1,v_relative.^2)/Ttotal;
    waitbar(i/nc,hw1,['第一种情况计算完成：',num2str(i/nc*100),'%'])
end
%%
figure
plot(C_constant_damp,pl,'r','Linewidth',2);
xlabel('定常阻尼系数');
ylabel('平均输出功率大小');
title('确定最优定常阻尼系数');
MaxP1=max(P1);
c=C_constant_damp(P1==MaxP1);
nc=100;
nmi=10;
[C0,MI]=meshgrid(linspace(0,100000,nc),linspace(0,1,nmi));
P2=zeros(size(C0));
hw2=waitbar(0,'第二种情况开始计算...');
pause(1)

for i=1:nmi
    for j=1:nc
        c0=C0(i,j);
        mi=MI(i,j);
        [tt2,yy2]=ode45(@(t,y)myode2(t,y,c0,mi),tspan,y0);
        v_relative=yy2(:,4)-yy2(:,4);
        c=c0*(abs(v_relative)).^mi;
        P2(i,j)=1/2*trapz(tt2,c.*v_relative.^2)/Ttotal;
        waitbar(((i-1)*nc+j)/(nc*nmi),hw2,['第二种情况计算完成：',num2str((i-1)*nc+j)/(nc*nmi)*100,'%'])
    end
end

%%
figure
surf(C0,MI,P2);
shading interp
camlight
lighting phong
xlabel('比例系数');
ylabel('幂指数');
zlabel('平均输出功率大小');
title('确定最优的非定常阻尼系数')
MaxP2=max(max(P2));
c0=C0(P2==MaxP2);
mi=MI(P2==MaxP2);

function dy=myode1(t,y,c)
global m1 m2 mf cx k omega f rho g r
dy=zeros(4,1);
dy(1)=y(2);
dy(2)=-cx/(m1+mf)*y(2)-pi*r^2*rho*g/(m1+mf)*y(1)-c/(m1+mf)*(y(2)-y(4))-k/(m1+mf)*(y(1)-y(3))+f/(m1+mf)*cos(omega*t);
dy(3)=y(4);
dy(4)=-c/m2*(y(4)-y(2))-k/m2*(y(3)-y(1));
end

function dy=myode2(t,y,c0,mi)
global m1 m2 mf cx k omega f rho g r
c=c0*(abs(y(2)-y(4)))^mi;
dy=zeros(4,1);
dy(1)=y(2);
dy(2)=-cx/(m1+mf)*y(2)-pi*r^2*rho*g/(m1+mf)*y(1)-c/(m1+mf)*(y(2)-y(4))-k/(m1+mf)*(y(1)-y(3))+f/(m1+mf)*cos(omega*t);
dy(3)=y(4);
dy(4)=-c/m2*(y(4)-y(2))-k/m2*(y(3)-y(1));
end





