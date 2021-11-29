 clear all;clf;cla;
 pi=3.14159;%Բ����
 rhoa=1.181;%���ڿ����ܶ�,��Ӧ26��
 u0=0.1;%���������ٶ�
 b0=0.2;%��ڿ��
 nu=1.48e-5;%�˶�ճ��
 g=-9.81;%�������ٶ�
 lambda=1;%ʵ�龭�����
 a=0.04;%����ϵ��
 rho00=1.205;%���ڿ����ܶ�
 theta=0;%���������Ƕ�
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

t0=0;%�Ա�����ֵ
%����x0��dx0�е���Ŀ�����ĳ�ֵ����ʵд��û�и���������д
x0=[u(i) b0 rho00 theta 0 1];%�����ֵx0�Ķ������ˣ����Ա��������
fix_x0=[1 1 1 1 1 1 1];%��ֵx0����Ŀ�����ģ���Ӧλ����1������Ϊ0,������x0�������ˣ���ȫΪ1
dx0=[0 0 0 0 0 0]';%�����г�ֵdx0һ����û�и�������ôȫ������д 
fix_dx0=zeros(6,1);%��ֵdx0����Ŀ�����ģ���Ӧλ����1������Ϊ0 ������dx0һ��û�и�������ȫ��Ϊ0
[x02,dx02]=decic(odefun,t0,x0,fix_x0,dx0,fix_dx0);
[t,x,dx]=ode15i(odefun,[0 10],x02,dx02);%x02��dx02����decic�������
 %--------------------------------------------------------------------------
 %���ò�ֵ�������ݶ��룬�����Ȳ�ͬ��һά����ͳһת��Ϊ200���ȵ�һά����
 L=length(x(:,2));%����x����������
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
xb=xb/D;%�����ٻ�
yb=yb/D;%�����ٻ�
b=b/D;%�����ٻ�
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




