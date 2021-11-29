 clear all;clf;cla;
 pi=3.14159;%圆周率
 rhoa=1.181;%室内空气密度,对应26℃
 u0=0.1;%射流出口速度
 b0=0.2;%风口宽度
 nu=1.48e-5;%运动粘度
 g=-9.81;%重力加速度
 lambda=1;%实验经验参数
 a=0.04;%卷吸系数
 rho00=1.205;%出口空气密度
 theta=0;%出口射流角度
 D=1;
 for i=1:5
    u(i)=i/10;

 Frd0(i)=u(i)/(sqrt(((rhoa-rho00)/rhoa)*g*2*b0));
 
 odefun=@(t,x,dx)[x(2)*dx(1)+x(1)*dx(2)-2*a*x(1)/sqrt(pi);
                  x(2)*(x(3)-rhoa)*dx(1)+x(1)*(x(3)-rhoa)*dx(2)+x(1)*x(2)*dx(3);
                  (x(1)*x(2)*x(3)*sqrt(pi)*cos(x(4))/sqrt(2))*dx(1)+(x(1)^2*x(3)*sqrt(pi)*cos(x(4))/(2*sqrt(2)))*dx(2)+(x(1)^2*x(2)*sqrt(pi)*cos(x(4))/(2*sqrt(2)))*dx(3)-(x(1)^2*x(2)*x(3)*sqrt(pi)*sin(x(4))/(2*sqrt(2)))*dx(4)+((0.2*x(3)*u(i)^2/2)*cos(x(4))/((x(6)/b0)*(u(i)*b0/nu)^1/12));
                  (x(1)*x(2)*x(3)*sqrt(pi)*sin(x(4))/sqrt(2))*dx(1)+(x(1)^2*x(3)*sqrt(pi)*sin(x(4))/(2*sqrt(2)))*dx(2)+(x(1)^2*x(2)*sqrt(pi)*sin(x(4))/(2*sqrt(2)))*dx(3)+(x(1)^2*x(2)*x(3)*sqrt(pi)*cos(x(4))/(2*sqrt(2)))*dx(4)+((0.2*x(3)*u(i)^2/2)*sin(x(4))/((x(6)/b0)*(u(i)*b0/nu)^1/12))-g*lambda*x(2)*sqrt(pi)*(x(3)-rhoa)/2;
                  dx(5)-cos(x(4));
                  dx(6)-sin(x(4))]; 

t0=0;%自变量初值
%对于x0和dx0中的题目给出的初值，如实写，没有给出的任意写
x0=[u(i) b0 rho00 theta 0 1];%本题初值x0的都给出了，所以必须是这个
fix_x0=[1 1 1 1 1 1 1];%初值x0中题目给出的，对应位置填1，否则为0,本题中x0都给出了，故全为1
dx0=[0 0 0 0 0 0]';%本题中初值dx0一个都没有给出，那么全部任意写 
fix_dx0=zeros(6,1);%初值dx0中题目给出的，对应位置填1，否则为0 本题中dx0一个没有给出，故全部为0
[x02,dx02]=decic(odefun,t0,x0,fix_x0,dx0,fix_dx0);
[t,x,dx]=ode15i(odefun,[0 10],x02,dx02);%x02和dx02是由decic输出参数
 %--------------------------------------------------------------------------
 %采用差值法将数据对齐，将长度不同的一维数组统一转化为200长度的一维数组
 L=length(x(:,2));%测量x变量的行数
 b(1:L,1)=x(:,2);
 b(200,1)=b(L,1);
 xb(1:L,1)=x(:,5);
 xb(200,1)=xb(L,1);
 yb(1:L,1)=x(:,6);
 yb(200,1)=yb(L,1);
 for j=1:201-L
     b(L-1+j,1)=b(L-1,1)+(b(200,1)-b(L-1,1))*j/(201-L);
     xb(L-1+j,1)=xb(L-1,1)+(xb(200,1)-xb(L-1,1))*j/(201-L);
     yb(L-1+j,1)=yb(L-1,1)+(yb(200,1)-yb(L-1,1))*j/(201-L);
 end
xb=xb/D;%无量纲化
yb=yb/D;%无量纲化
b=b/D;%无量纲化
for k=1:199
     if k==1 
         ds(1)=sqrt((xb(2)-xb(1))^2+(yb(2)-yb(1))^2);
         S(1)=ds(1);
     end
     if k~=1
         ds(k)=sqrt((xb(k+1)-xb(k))^2+(yb(k+1)-yb(k))^2);
         S(k)=S(k-1)+ds(k);
     end

end
 S(200)=S(199);
%--------------------------------------------------------------------------
% if i==1
%     plot(b(:,1),S(1,:),'black','LineWidth',1);
% end
% hold on
% if i==2
%     plot(b(:,1),S(1,:),'black','LineWidth',1);
% end
% hold on
% if i==3
%     plot(b(:,1),S(1,:),'black','LineWidth',1);
% end
% hold on
% if i==4
%     plot(b(:,1),S(1,:),'black','LineWidth',1);
% end
% hold on
% if i==5
%      plot(b(:,1),S(1,:),'black','LineWidth',1);
% end
% hold on
if i==1
    plot(x(:,2),x(:,6),'black','LineWidth',1);
end
hold on
if i==2
   plot(x(:,2),x(:,6),'black','LineWidth',1);
end
hold on
if i==3
    plot(x(:,2),x(:,6),'black','LineWidth',1);
end
hold on
if i==4
    plot(x(:,2),x(:,6),'black','LineWidth',1);
end
hold on
if i==5
   plot(x(:,2),x(:,6),'black','LineWidth',1);
end
hold on
 end
 
grid minor
set(gca,'ygrid','on','xgrid','on','GridLineStyle','-','GridAlpha',0.8,'MinorGridLineStyle','-','MinorGridAlpha',0.2)
set(gca,'xaxislocation','bottom');
set(gca,'Fontname','Times New Roman');
ylabel('S')
xlabel('b')
% set(gca,'YLim',[0 1]);




