function J1=compute_J()
m1=4866;
ha=2/3*0.8;
hb=2.3;
hc=3.8;
Sa=pi*sqrt(1+0.8^2);
Sb=6*pi;
Sc=pi;
Sabc=Sa+Sb+Sc;
zc=Sa/Sabc*ha+Sb/Sabc*hb+Sc/Sabc*hc;
rho=m1/Sabc; % 密度
Ja=rho*sqrt(1+0.8^2)*(1/2+0.8^2)/4*2*pi+rho*sqrt(1+0.8^2)*pi*zc^2-2*rho*sqrt(1+0.8^2)*zc*0.8/3*2*pi; % 圆锥部分 
Jb=rho*(pi*3+2*pi*zc^2*3+2*pi/3*54.36-2*pi*zc*13.8); % 圆柱部分 
Jc=rho*(3.8-zc)^2*pi+rho*2*pi/2/4; % 圆盖部分 
J1=Ja+Jb+Jc;  % 总的对质心转动惯量
end