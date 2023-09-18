clc,clear,close all
%%
%表3、图1
%基准参数设定
theta=0.5;sigma=0.1;e=0.8;belta=0.33;gamma=0.15;omega=1;epsilon=0.4;alpha=0.25;delta=1;
n=1.14;g=2.5;mu=0.8;
l=0.4;
%基准参数校准
u=epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega);
syms R;
F=@(mu,theta,R,alpha,delta,l,e,belta,n,g,u)@(R)(1-R*mu*theta/(R+belta*R+(1-e)/e*mu*omega))*(1+u)-((1-R*mu*theta/(R+belta*R+(1-e)/e*mu*omega))+R*(omega-theta)*mu);
R=fzero(F(mu,theta,R,alpha,delta,l,e,belta,n,g,u),rand);
syms v
F2=@(R,omega,theta,mu,e,belta,n,v)@(v)R*(omega-theta)*mu/(((l-R*mu*theta/((belta*R*mu+(1-e)/e*mu*omega)+R))*(belta*R*mu+(1-e)/e*mu*omega+R)+R*(omega-theta)*mu))*(1+u)/v-n;
v=fzero(F2(R,omega,theta,mu,e,belta,n,v),rand);
B=log(g)*(1+sigma)-theta/(1+sigma)*log(1/e*theta*v/(omega-theta)*(1+l*delta/R))-sigma/(1+sigma)*log(delta*(1-v*n+u));
%基准数值模拟
%e=0.8;
l=0.4:0.01:0.6;
%家庭隔代照料时间
FU=epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega);
for i=1:length(l)
FU(i)=epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega);
end
U = [FU];
%青年人就业机会
FN=R*mu*(omega-theta)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+R*mu*(omega-theta))*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega))/v;
for i=1:length(l)
FN(i)=R*mu*(omega-theta)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+R*mu*(omega-theta))*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega))/v;
end
N = [FN];
%青年人就业意愿
FVN=(1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega));
for i=1:length(l)
FVN(i)=(1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega));
end
VN = [FVN];
%青年人就业能力
FG=1/(1+sigma)*B+log(1/e*theta*v/(omega-theta)*(1+l*delta/R))*theta/(1+sigma)+log(delta*((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega))))*sigma/(1+sigma);
for i=1:length(l)
FG(i)=1/(1+sigma)*B+log(1/e*theta*v/(omega-theta)*(1+l(i)*delta/R))*theta/(1+sigma)+log(delta*((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega))))*sigma/(1+sigma);
end
G = exp(FG);

%e=0.5
e=0.5;
l=0.4:0.01:0.6;

FU1=epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega);
for i=1:length(l)
FU1(i)=epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega);
end
U1 = [FU1];

FN1=R*mu*(omega-theta)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+R*mu*(omega-theta))*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega))/v+1;
for i=1:length(l)
FN1(i)=R*mu*(omega-theta)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+R*mu*(omega-theta))*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega))/v+1;
end
N1 = [FN1];

FVN1=(1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega));
for i=1:length(l)
FVN1(i)=(1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega));
end
VN1 = [FVN1];

FG1=1/(1+sigma)*B+log(1/e*theta*v/(omega-theta)*(1+l*delta/R))*theta/(1+sigma)+log(delta*((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega))))*sigma/(1+sigma);
for i=1:length(l)
FG1(i)=1/(1+sigma)*B+log(1/e*theta*v/(omega-theta)*(1+l(i)*delta/R))*theta/(1+sigma)+log(delta*((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega))))*sigma/(1+sigma);
end

G1 = exp(FG1);
%e=0.2
e=0.2;
l=0.4:0.01:0.6;

FU2=epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega);
for i=1:length(l)
FU2(i)=epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega);
end
U2 = [FU2];

FN2=R*mu*(omega-theta)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+R*mu*(omega-theta))*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega))/v+1;
for i=1:length(l)
FN2(i)=R*mu*(omega-theta)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+R*mu*(omega-theta))*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega))/v+1;
end
N2 = [FN2];

FVN2=(1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega));
for i=1:length(l)
FVN2(i)=(1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega));
end
VN2 = [FVN2];

FG2=1/(1+sigma)*B+log(1/e*theta*v/(omega-theta)*(1+l*delta/R))*theta/(1+sigma)+log(delta*((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega))))*sigma/(1+sigma);
for i=1:length(l)
FG2(i)=1/(1+sigma)*B+log(1/e*theta*v/(omega-theta)*(1+l(i)*delta/R))*theta/(1+sigma)+log(delta*((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega))))*sigma/(1+sigma);
end
G2 = exp(FG2);
%作图
clf;
subplot(2,2,1),plot(l,U,'-',l,U1,'--',l,U2,'-.');
subplot(2,2,2),plot(l,N,'-',l,N1,'--',l,N2,'-.');
subplot(2,2,3),plot(l,VN,'-',l,VN1,'--',l,VN2,'-.');
subplot(2,2,4),plot(l,G,'-',l,G1,'--',l,G2,'-.');

subplot(2,2,1);
set(gca,'linewidth',1.5); 
set(gca,'xtick',[0.4,0.44,0.48,0.52,0.56,0.6]);
set(gca,'xticklabel',{'60','61','62','63','64','65'});
set(gca,'fontsize',12);
set (gca,'FontSize',12)
xlabel('退休年龄')
ylabel('隔代照料程度')

subplot(2,2,2);
set(gca,'linewidth',1.5); 
set(gca,'xtick',[0.4,0.44,0.48,0.52,0.56,0.6]);
set(gca,'xticklabel',{'60','61','62','63','64','65'});
set(gca,'fontsize',12);
set (gca,'FontSize',12)
xlabel('退休年龄')
ylabel('青年人就业机会')

subplot(2,2,3);
set(gca,'linewidth',1.5); 
set(gca,'xtick',[0.4,0.44,0.48,0.52,0.56,0.6]);
set(gca,'xticklabel',{'60','61','62','63','64','65'});
set(gca,'fontsize',12);
set (gca,'FontSize',12)
xlabel('退休年龄')
ylabel('青年人就业意愿')

subplot(2,2,4);
set(gca,'linewidth',1.5); 
set(gca,'xtick',[0.4,0.44,0.48,0.52,0.56,0.6]);
set(gca,'xticklabel',{'60','61','62','63','64','65'});
set(gca,'fontsize',12);
set (gca,'FontSize',12)
xlabel('退休年龄')
ylabel('青年人就业能力')
%%
clc,clear,close all
% 图2
%基准参数设定
theta=0.5;sigma=0.1;e=0.8;belta=0.33;gamma=0.15;omega=1;epsilon=0.4;alpha=0.25;delta=1;
n=1.14;g=2.5;mu=0.8;
l=0.4;t=0.24;
%基准参数校准
syms R;
F=@(epsilon,theta,belta,alpha,R,t,delta,n,g)@(R)((1-t)*R-t*n-(1+mu*theta)/(1+belta+epsilon*theta+mu*theta)*((1-t)*R+delta*l))-alpha/(1-alpha)*(n*g+delta*l);
R=fzero(F(epsilon,theta,belta,alpha,R,t,delta,n,g),rand);
syms D;
D=t*g*n*(1+epsilon*theta+belta+mu*theta)/(delta*((1-t)*R+l*delta))+epsilon*omega;
syms u;
u=D/(D+gamma)*(1-l)-gamma/(D+gamma);
syms v;
v=mu*(omega-theta)/(1+epsilon*theta+belta+mu*omega)*(1+u)/n;
syms B;
B=g^(1+delta)*(1/e*v/(omega-theta)*((D-epsilon)*delta/mu+theta)*(1-t+delta*l/R))^(-theta)*(delta*((1+epsilon*theta+belta+mu*theta)/(1+epsilon*theta+belta+mu*omega)*(1+u)))^(-delta);
l=0.4:0.01:0.6;%延迟退休到65；
%基准数值模拟；
FU=D/(D+gamma)*(1-l)-gamma/(D+gamma);
for i=1:length(l)
FU(i)=D/(D+gamma)*(1-l(i))-gamma/(D+gamma);
end
U = [FU];

FN=1/v*mu*(omega-theta)/(1+epsilon*theta+belta+mu*theta)*(1+D/(D+gamma)*(1-l)-gamma/(D+gamma));
for i=1:length(l)
FN(i)=1/v*mu*(omega-theta)/(1+epsilon*theta+belta+mu*theta)*(1+D/(D+gamma)*(1-l(i))-gamma/(D+gamma));
end
N = [FN];

FVN=(1+epsilon*theta+belta+mu*theta)/(1+epsilon*theta+belta+mu*omega)*(1+D/(D+gamma)*(1-l)-gamma/(D+gamma));
for i=1:length(l)
FVN(i)=(1+epsilon*theta+belta+mu*theta)/(1+epsilon*theta+belta+mu*omega)*(1+D/(D+gamma)*(1-l(i))-gamma/(D+gamma));

end
VN = [FVN];

FG =log(B)*(1/(1+sigma))+log(1/e*v/(omega-theta)*((D-epsilon)*delta/mu+theta)*(1-t+delta*l/R))*(theta/(1+sigma))+log(delta*((1+epsilon*theta+belta+mu*theta)/(1+epsilon*theta+belta+mu*omega)*(1+D/(D+gamma)*(1-l)-gamma/(D+gamma))))*(sigma/(1+sigma));
for i=1:length(l)
FG(i)=log(B)*(1/(1+sigma))+log(1/e*v/(omega-theta)*((D-epsilon)*delta/mu+theta)*(1-t+delta*l(i)/R))*(theta/(1+sigma))+log(delta*((1+epsilon*theta+belta+mu*theta)/(1+epsilon*theta+belta+mu*omega)*(1+D/(D+gamma)*(1-l(i))-gamma/(D+gamma))))*(sigma/(1+sigma));
end
H= [FG];
G=exp(H);
%现收现付制参数设定
t=0.16;l=0.4;
%现收现付制参数校准
syms R;
F=@(e,epsilon,theta,belta,alpha,R,t,delta,n,g)@(R)-(1/(1+e/(1-e)*epsilon*theta/belta)+1/((1+epsilon)*belta))*alpha/(1+alpha)+((1-t)*(R)/(1+e/(1-e)*epsilon*theta/belta)+t/delta*n/((1+epsilon*theta)*belta))/(n*g+delta);
R=fzero(F(e,epsilon,theta,belta,alpha,R,t,delta,n,g),rand);
syms D;
D=t*g*n*(1+epsilon*theta+belta+mu*theta)/(delta*((1-t)*R+l*delta))+epsilon*omega;
syms u;
u=D/(D+gamma)*(1-l)-gamma/(D+gamma);
syms v;
v=mu*(omega-theta)/(1+epsilon*theta+belta+mu*omega)*(1+u)/n;
syms B;
B=g^(1+delta)*(1/e*v/(omega-theta)*((D-epsilon)*delta/mu+theta)*(1-t+delta*l/R))^(-theta)*(delta*((1+epsilon*theta+belta+mu*theta)/(1+epsilon*theta+belta+mu*omega)*(1+u)))^(-delta);
%现收现付制数值模拟
l=0.4:0.01:0.6;%延迟退休到65；
FU1=D/(D+gamma)*(1-l)-gamma/(D+gamma);
for i=1:length(l)
FU1(i)=D/(D+gamma)*(1-l(i))-gamma/(D+gamma);
end
U1 = [FU1];

FN1=1/v*mu*(omega-theta)/(1+epsilon*theta+belta+mu*theta)*(1+D/(D+gamma)*(1-l)-gamma/(D+gamma));
for i=1:length(l)
FN1(i)=1/v*mu*(omega-theta)/(1+epsilon*theta+belta+mu*theta)*(1+D/(D+gamma)*(1-l(i))-gamma/(D+gamma));
end
N1 = [FN1];

FVN1=(1+epsilon*theta+belta+mu*theta)/(1+epsilon*theta+belta+mu*omega)*(1+D/(D+gamma)*(1-l)-gamma/(D+gamma));
for i=1:length(l)
FVN1(i)=(1+epsilon*theta+belta+mu*theta)/(1+epsilon*theta+belta+mu*omega)*(1+D/(D+gamma)*(1-l(i))-gamma/(D+gamma));
end
VN1 = [FVN1];

FG1 =log(B)*(1/(1+sigma))+log(1/e*v/(omega-theta)*((D-epsilon)*delta/mu+theta)*(1-t+delta*l/R))*(theta/(1+sigma))+log(delta*((1+epsilon*theta+belta+mu*theta)/(1+epsilon*theta+belta+mu*omega)*(1+D/(D+gamma)*(1-l)-gamma/(D+gamma))))*(sigma/(1+sigma));
for i=1:length(l)
FG1(i)=log(B)*(1/(1+sigma))+log(1/e*v/(omega-theta)*((D-epsilon)*delta/mu+theta)*(1-t+delta*l(i)/R))*(theta/(1+sigma))+log(delta*((1+epsilon*theta+belta+mu*theta)/(1+epsilon*theta+belta+mu*omega)*(1+D/(D+gamma)*(1-l(i))-gamma/(D+gamma))))*(sigma/(1+sigma));
end
H1= [FG1];
G1=exp(H1);
%作图
clf;
subplot(2,2,1),plot(l,U,'-',l,U1,'--');
subplot(2,2,2),plot(l,N,'-',l,N1,'--');
subplot(2,2,3),plot(l,VN,'-',l,VN1,'--');
subplot(2,2,4),plot(l,G,'-',l,G1,'--');

subplot(2,2,1);
set(gca,'linewidth',1.5); 
set(gca,'xtick',[0.4,0.44,0.48,0.52,0.56,0.6]);
set(gca,'xticklabel',{'60','61','62','63','64','65'});
set(gca,'fontsize',12);
set (gca,'FontSize',12)
xlabel('退休年龄')
ylabel('隔代照料程度')

subplot(2,2,2);
set(gca,'linewidth',1.5); 
set(gca,'xtick',[0.4,0.44,0.48,0.52,0.56,0.6]);
set(gca,'xticklabel',{'60','61','62','63','64','65'});
set(gca,'fontsize',12);
set (gca,'FontSize',12)
xlabel('退休年龄')
ylabel('青年人就业机会')

subplot(2,2,3);
set(gca,'linewidth',1.5); 
set(gca,'xtick',[0.4,0.44,0.48,0.52,0.56,0.6]);
set(gca,'xticklabel',{'60','61','62','63','64','65'});
set(gca,'fontsize',12);
set (gca,'FontSize',12)
xlabel('退休年龄')
ylabel('青年人就业意愿')

subplot(2,2,4);
set(gca,'linewidth',1.5); 
set(gca,'xtick',[0.4,0.44,0.48,0.52,0.56,0.6]);
set(gca,'xticklabel',{'60','61','62','63','64','65'});
set(gca,'fontsize',12);
set (gca,'FontSize',12)
xlabel('退休年龄')
ylabel('青年人就业能力')

%%
clc,clear,close all
% 图3
%图3
clc,clear,close all
%参数设定
theta=0.5;sigma=0.1;e=0.8;belta=0.33;gamma=0.15;omega=1;epsilon=0.4;alpha=0.25;delta=1;
n=1.14;g=2.5;mu=0.8;
l=0.4;t=0.24;
%现收现付参数校准
syms R;
F=@(epsilon,theta,belta,alpha,R,t,delta,n,g)@(R)((1-t)*R-t*n-(1+mu*theta)/(1+belta+epsilon*theta+mu*theta)*((1-t)*R+delta*l))-alpha/(1-alpha)*(n*g+delta*l);
R=fzero(F(epsilon,theta,belta,alpha,R,t,delta,n,g),rand);
syms D;
D=t*g*n*(1+epsilon*theta+belta+mu*theta)/(delta*((1-t)*R+l*delta))+epsilon*omega;
syms u;
u=D/(D+gamma)*(1-l)-gamma/(D+gamma);
syms v;
v=mu*(omega-theta)/(1+epsilon*theta+belta+mu*omega)*(1+u)/n;
syms B;
B=g^(1+delta)*(1/e*v/(omega-theta)*((D-epsilon)*delta/mu+theta)*(1-t+delta*l/R))^(-theta)*(delta*((1+epsilon*theta+belta+mu*theta)/(1+epsilon*theta+belta+mu*omega)*(1+u)))^(-delta);
%现收现付制数值模拟；
l=0.4:0.01:0.6;%延迟退休到65；
FU=D/(D+gamma)*(1-l)-gamma/(D+gamma);
for i=1:length(l)
FU(i)=D/(D+gamma)*(1-l(i))-gamma/(D+gamma);
end
U = [FU];

FN=1/v*mu*(omega-theta)/(1+epsilon*theta+belta+mu*theta)*(1+D/(D+gamma)*(1-l)-gamma/(D+gamma));
for i=1:length(l)
FN(i)=1/v*mu*(omega-theta)/(1+epsilon*theta+belta+mu*theta)*(1+D/(D+gamma)*(1-l(i))-gamma/(D+gamma));
end
N = [FN];

FVN=(1+epsilon*theta+belta+mu*theta)/(1+epsilon*theta+belta+mu*omega)*(1+D/(D+gamma)*(1-l)-gamma/(D+gamma));
for i=1:length(l)
FVN(i)=(1+epsilon*theta+belta+mu*theta)/(1+epsilon*theta+belta+mu*omega)*(1+D/(D+gamma)*(1-l(i))-gamma/(D+gamma));

end
VN = [FVN];

FG =log(B)*(1/(1+sigma))+log(1/e*v/(omega-theta)*((D-epsilon)*delta/mu+theta)*(1-t+delta*l/R))*(theta/(1+sigma))+log(delta*((1+epsilon*theta+belta+mu*theta)/(1+epsilon*theta+belta+mu*omega)*(1+D/(D+gamma)*(1-l)-gamma/(D+gamma))))*(sigma/(1+sigma));
for i=1:length(l)
FG(i)=log(B)*(1/(1+sigma))+log(1/e*v/(omega-theta)*((D-epsilon)*delta/mu+theta)*(1-t+delta*l(i)/R))*(theta/(1+sigma))+log(delta*((1+epsilon*theta+belta+mu*theta)/(1+epsilon*theta+belta+mu*omega)*(1+D/(D+gamma)*(1-l(i))-gamma/(D+gamma))))*(sigma/(1+sigma));
end
H= [FG];
G=exp(H);
%基金值参数设定
theta=0.4;sigma=0.1;e=0.1;belta=0.33;gamma=0.15;omega=1;epsilon=0.4;alpha=0.25;delta=1;
n=1.14;g=2.5;mu=0.8;
l=0.4;t=0.24;

%基金值参数校准
A=1+(1-e)/e*mu*theta;
syms v;
v=1/n*mu*(omega-theta)/(1+A*belta+mu*omega)*epsilon/(gamma+epsilon)*(2-l);
syms R;
F=@(alpha,n,g,delta,mu,theta,A,belta,l,t)@(R)alpha/(1-alpha)-1/(n*g+delta)*(R-(1+mu*theta)/(1+mu*theta+A*belta)*(R+delta*(1+t)*l));
R=fzero(F(alpha,n,g,delta,mu,theta,A,belta,l,t),rand);
syms B;
B=g^(1+delta)*(1/e*v*theta/(omega-theta)*(1+delta/R*(1+t)*l))^(-theta)*(delta*(1+A*belta+mu*theta)/(1+A*belta+mu*omega)*epsilon/(gamma+epsilon)*(2-l))^(-delta);
%基金制数值模拟
l=0.4:0.01:0.6;%延迟退休到65；
FUj=epsilon/(gamma+epsilon)*(1-l)-gamma/(gamma+epsilon);
for i=1:length(l)
FUj(i)=epsilon/(gamma+epsilon)*(1-l(i))-gamma/(gamma+epsilon);
end
Uj = [FUj];

FNj=mu*(omega-theta)/(1+A*belta+mu*omega)*1/v*epsilon/(gamma+epsilon)*(2-l);
for i=1:length(l)
FNj(i)=mu*(omega-theta)/(1+A*belta+mu*omega)*1/v*epsilon/(gamma+epsilon)*(2-l(i));
end
Nj = [FNj];


FVNj = (1+A*belta+mu*theta)/(1+A*belta+mu*omega)*(epsilon/(gamma+epsilon))*(2-l);
for i=1:length(l)
FVNj(i)= (1+A*belta+mu*theta)/(1+A*belta+mu*omega)*(epsilon/(gamma+epsilon))*(2-l(i));
end
VNj = [FVNj];

FGj =1/(1+sigma)*log(B)+theta/(1+sigma)*log(1/e*v*theta/(omega-theta)*(1+delta/R*(1+t)*l))+sigma/(1+sigma)*log(delta*(1+A*belta+mu*theta)/(1+A*belta+mu*omega)*(epsilon/(gamma+epsilon))*(2-l));
for i=1:length(l)
FGj(i)=1/(1+sigma)*log(B)+theta/(1+sigma)*log(1/e*v*theta/(omega-theta)*(1+delta/R*(1+t)*l(i)))+sigma/(1+sigma)*log(delta*(1+A*belta+mu*theta)/(1+A*belta+mu*omega)*(epsilon/(gamma+epsilon))*(2-l(i)));
end
Hj= [FGj];
Gj=exp(FGj);
%作图
clf;
subplot(2,2,1),plot(l,U,'-',l,Uj,'--');
subplot(2,2,2),plot(l,N,'-',l,Nj,'--');
subplot(2,2,3),plot(l,VN,'-',l,VNj,'--');
subplot(2,2,4),plot(l,G,'-',l,Gj,'--');

subplot(2,2,1);
set(gca,'linewidth',1.5); 
set(gca,'xtick',[0.4,0.44,0.48,0.52,0.56,0.6]);
set(gca,'xticklabel',{'60','61','62','63','64','65'});
set(gca,'fontsize',12);
set (gca,'FontSize',12)
xlabel('退休年龄')
ylabel('隔代照料程度')

subplot(2,2,2);
set(gca,'linewidth',1.5); 
set(gca,'xtick',[0.4,0.44,0.48,0.52,0.56,0.6]);
set(gca,'xticklabel',{'60','61','62','63','64','65'});
set(gca,'fontsize',12);
set (gca,'FontSize',12)
xlabel('退休年龄')
ylabel('青年人就业机会')

subplot(2,2,3);
set(gca,'linewidth',1.5); 
set(gca,'xtick',[0.4,0.44,0.48,0.52,0.56,0.6]);
set(gca,'xticklabel',{'60','61','62','63','64','65'});
set(gca,'fontsize',12);
set (gca,'FontSize',12)
xlabel('退休年龄')
ylabel('青年人就业意愿')

subplot(2,2,4);
set(gca,'linewidth',1.5); 
set(gca,'xtick',[0.4,0.44,0.48,0.52,0.56,0.6]);
set(gca,'xticklabel',{'60','61','62','63','64','65'});
set(gca,'fontsize',12);
set (gca,'FontSize',12)
xlabel('退休年龄')
ylabel('青年人就业能力')
%%
%附表6
clc,clear,close all
%参数设定
theta=0.5;sigma=0.1;e=0.8;belta=0.33;gamma=0.15;omega=1;epsilon=0.4;alpha=0.25;delta=1;
n=1.14;g=2.5;mu=0.8;
l=0.4;
u=epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega);
syms R;
F=@(mu,theta,R,alpha,delta,l,e,belta,n,g,u)@(R)(1-R*mu*theta/(R+belta*R+(1-e)/e*mu*omega))*(1+u)-((1-R*mu*theta/(R+belta*R+(1-e)/e*mu*omega))+R*(omega-theta)*mu);
R=fzero(F(mu,theta,R,alpha,delta,l,e,belta,n,g,u),rand);
syms v
F2=@(R,omega,theta,mu,e,belta,n,v)@(v)R*(omega-theta)*mu/(((l-R*mu*theta/((belta*R*mu+(1-e)/e*mu*omega)+R))*(belta*R*mu+(1-e)/e*mu*omega+R)+R*(omega-theta)*mu))*(1+u)/v-n;
v=fzero(F2(R,omega,theta,mu,e,belta,n,v),rand);
B=log(g)*(1+sigma)-theta/(1+sigma)*log(1/e*theta*v/(omega-theta)*(1+l*delta/R))-sigma/(1+sigma)*log(delta*(1-v*n+u));

% gamma=0.15;sigma=0.1
l=0.4:0.01:0.6;e=0.8;
FU=epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega);
for i=1:length(l)
FU(i)=epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega);
end
U = [FU];

FN=R*mu*(omega-theta)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+R*mu*(omega-theta))*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega))/v;
for i=1:length(l)
FN(i)=R*mu*(omega-theta)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+R*mu*(omega-theta))*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega))/v;
end
N = [FN];

FVN=(1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega));
for i=1:length(l)
FVN(i)=(1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega));
end
VN = [FVN];

FG=1/(1+sigma)*B+log(1/e*theta*v/(omega-theta)*(1+l*delta/R))*theta/(1+sigma)+log(delta*((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega))))*sigma/(1+sigma);
for i=1:length(l)
FG(i)=1/(1+sigma)*B+log(1/e*theta*v/(omega-theta)*(1+l(i)*delta/R))*theta/(1+sigma)+log(delta*((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega))))*sigma/(1+sigma);
end
G = exp(FG);
%sigma=0.04;gamma=0.15;
sigma=0.04;gamma=0.15;
FU1=epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega);
for i=1:length(l)
FU1(i)=epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega);
end
U1 = [FU1];

FN1=R*mu*(omega-theta)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+R*mu*(omega-theta))*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega))/v;
for i=1:length(l)
FN1(i)=R*mu*(omega-theta)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+R*mu*(omega-theta))*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega))/v;
end
N1 = [FN1];

FVN1=(1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega));
for i=1:length(l)
FVN1(i)=(1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega));
end
VN1 = [FVN1];


FG1=1/(1+sigma)*B+log(1/e*theta*v/(omega-theta)*(1+l*delta/R))*theta/(1+sigma)+log(delta*((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega))))*sigma/(1+sigma);
for i=1:length(l)
FG1(i)=1/(1+sigma)*B+log(1/e*theta*v/(omega-theta)*(1+l(i)*delta/R))*theta/(1+sigma)+log(delta*((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega))))*sigma/(1+sigma);
end

G1 = exp(FG1);


%sigma=0.16;gamma=0.15;
sigma=0.16;gamma=0.15;
FU2=epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega);
for i=1:length(l)
FU2(i)=epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega);
end
U2 = [FU2];


FN2=R*mu*(omega-theta)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+R*mu*(omega-theta))*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega))/v;
for i=1:length(l)
FN2(i)=R*mu*(omega-theta)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+R*mu*(omega-theta))*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega))/v;
end
N2 = [FN2];

FVN2=(1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega));
for i=1:length(l)
FVN1(i)=(1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega));
end
VN2 = [FVN2];


FG2=1/(1+sigma)*B+log(1/e*theta*v/(omega-theta)*(1+l*delta/R))*theta/(1+sigma)+log(delta*((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega))))*sigma/(1+sigma);
for i=1:length(l)
FG2(i)=1/(1+sigma)*B+log(1/e*theta*v/(omega-theta)*(1+l(i)*delta/R))*theta/(1+sigma)+log(delta*((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega))))*sigma/(1+sigma);
end

G2 = exp(FG2);



%sigma=0.1;gamma=0.1;
sigma=0.1;gamma=0.1;
FU3=epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega);
for i=1:length(l)
FU3(i)=epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega);
end
U3 = [FU3];


FN3=R*mu*(omega-theta)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+R*mu*(omega-theta))*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega))/v;
for i=1:length(l)
FN3(i)=R*mu*(omega-theta)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+R*mu*(omega-theta))*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega))/v;
end
N3 = [FN3];

FVN3=(1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega));
for i=1:length(l)
FVN3(i)=(1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega));
end
VN3 = [FVN3];


FG3=1/(1+sigma)*B+log(1/e*theta*v/(omega-theta)*(1+l*delta/R))*theta/(1+sigma)+log(delta*((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega))))*sigma/(1+sigma);
for i=1:length(l)
FG3(i)=1/(1+sigma)*B+log(1/e*theta*v/(omega-theta)*(1+l(i)*delta/R))*theta/(1+sigma)+log(delta*((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega))))*sigma/(1+sigma);
end

G3 = exp(FG3);

%sigma=0.1;gamma=0.2;
sigma=0.1;gamma=0.2;
FU4=epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega);
for i=1:length(l)
FU4(i)=epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega);
end
U4 = [FU4];

FN4=R*mu*(omega-theta)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+R*mu*(omega-theta))*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega))/v;
for i=1:length(l)
FN4(i)=R*mu*(omega-theta)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+R*mu*(omega-theta))*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega))/v;
end
N4 = [FN4];

FVN4=(1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega));
for i=1:length(l)
FVN4(i)=(1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega));
end
VN4 = [FVN4];

FG4=1/(1+sigma)*B+log(1/e*theta*v/(omega-theta)*(1+l*delta/R))*theta/(1+sigma)+log(delta*((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega))))*sigma/(1+sigma);
for i=1:length(l)
FG4(i)=1/(1+sigma)*B+log(1/e*theta*v/(omega-theta)*(1+l(i)*delta/R))*theta/(1+sigma)+log(delta*((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega))))*sigma/(1+sigma);
end
G4 = exp(FG4);

%%
%附表7
clc,clear,close all
theta=0.4;sigma=0.1;e=0.1;belta=0.33;gamma=0.15;omega=1;epsilon=0.4;alpha=0.25;delta=1;
n=1.14;g=2.5;mu=0.8;
l=0.4;t=0.24;

%校准参数
A=1+(1-e)/e*mu*theta;
syms v;
v=1/n*mu*(omega-theta)/(1+A*belta+mu*omega)*epsilon/(gamma+epsilon)*(2-l);

syms R;
F=@(alpha,n,g,delta,mu,theta,A,belta,l,t)@(R)alpha/(1-alpha)-1/(n*g+delta)*(R-(1+mu*theta)/(1+mu*theta+A*belta)*(R+delta*(1+t)*l));
R=fzero(F(alpha,n,g,delta,mu,theta,A,belta,l,t),rand);

syms B;
B=g^(1+delta)*(1/e*v*theta/(omega-theta)*(1+delta/R*(1+t)*l))^(-theta)*(delta*(1+A*belta+mu*theta)/(1+A*belta+mu*omega)*epsilon/(gamma+epsilon)*(2-l))^(-delta);

l=0.4:0.01:0.6;%延迟退休到65；
FUj=epsilon/(gamma+epsilon)*(1-l)-gamma/(gamma+epsilon);
for i=1:length(l)
FUj(i)=epsilon/(gamma+epsilon)*(1-l(i))-gamma/(gamma+epsilon);
end
Uj = [FUj];

FNj=mu*(omega-theta)/(1+A*belta+mu*omega)*1/v*epsilon/(gamma+epsilon)*(2-l);
for i=1:length(l)
FNj(i)=mu*(omega-theta)/(1+A*belta+mu*omega)*1/v*epsilon/(gamma+epsilon)*(2-l(i));
end
Nj = [FNj];



FVNj = (1+A*belta+mu*theta)/(1+A*belta+mu*omega)*(epsilon/(gamma+epsilon))*(2-l);
for i=1:length(l)
FVNj(i)= (1+A*belta+mu*theta)/(1+A*belta+mu*omega)*(epsilon/(gamma+epsilon))*(2-l(i));
end
VNj = [FVNj];


%对青年人就业质量进行基准数值模拟；
FGj =1/(1+sigma)*log(B)+theta/(1+sigma)*log(1/e*v*theta/(omega-theta)*(1+delta/R*(1+t)*l))+sigma/(1+sigma)*log(delta*(1+A*belta+mu*theta)/(1+A*belta+mu*omega)*(epsilon/(gamma+epsilon))*(2-l));
for i=1:length(l)
FGj(i)=1/(1+sigma)*log(B)+theta/(1+sigma)*log(1/e*v*theta/(omega-theta)*(1+delta/R*(1+t)*l(i)))+sigma/(1+sigma)*log(delta*(1+A*belta+mu*theta)/(1+A*belta+mu*omega)*(epsilon/(gamma+epsilon))*(2-l(i)));
end

Hj= [FGj];

Gj=exp(FGj);


%sigma=0.04;gamma=0.15;
sigma=0.04;gamma=0.15;
l=0.4:0.01:0.6;%延迟退休到65；
FUj1=epsilon/(gamma+epsilon)*(1-l)-gamma/(gamma+epsilon);
for i=1:length(l)
FUj1(i)=epsilon/(gamma+epsilon)*(1-l(i))-gamma/(gamma+epsilon);
end
Uj1 = [FUj1];

FNj1=mu*(omega-theta)/(1+A*belta+mu*omega)*1/v*epsilon/(gamma+epsilon)*(2-l);
for i=1:length(l)
FNj1(i)=mu*(omega-theta)/(1+A*belta+mu*omega)*1/v*epsilon/(gamma+epsilon)*(2-l(i));
end
Nj1 = [FNj1];


FVNj1 = (1+A*belta+mu*theta)/(1+A*belta+mu*omega)*(epsilon/(gamma+epsilon))*(2-l);
for i=1:length(l)
FVNj1(i)= (1+A*belta+mu*theta)/(1+A*belta+mu*omega)*(epsilon/(gamma+epsilon))*(2-l(i));
end
VNj1 = [FVNj1];


FGj1 =1/(1+sigma)*log(B)+theta/(1+sigma)*log(1/e*v*theta/(omega-theta)*(1+delta/R*(1+t)*l))+sigma/(1+sigma)*log(delta*(1+A*belta+mu*theta)/(1+A*belta+mu*omega)*(epsilon/(gamma+epsilon))*(2-l));
for i=1:length(l)
FGj1(i)=1/(1+sigma)*log(B)+theta/(1+sigma)*log(1/e*v*theta/(omega-theta)*(1+delta/R*(1+t)*l(i)))+sigma/(1+sigma)*log(delta*(1+A*belta+mu*theta)/(1+A*belta+mu*omega)*(epsilon/(gamma+epsilon))*(2-l(i)));
end

Hj1= [FGj1];

Gj1=exp(FGj1);


%sigma=0.16;gamma=0.15;
sigma=0.16;gamma=0.15;
l=0.4:0.01:0.6;%延迟退休到65；
FUj2=epsilon/(gamma+epsilon)*(1-l)-gamma/(gamma+epsilon);
for i=1:length(l)
FUj2(i)=epsilon/(gamma+epsilon)*(1-l(i))-gamma/(gamma+epsilon);
end
Uj2 = [FUj2];

FNj2=mu*(omega-theta)/(1+A*belta+mu*omega)*1/v*epsilon/(gamma+epsilon)*(2-l);
for i=1:length(l)
FNj2(i)=mu*(omega-theta)/(1+A*belta+mu*omega)*1/v*epsilon/(gamma+epsilon)*(2-l(i));
end
Nj2 = [FNj2];


FVNj2 = (1+A*belta+mu*theta)/(1+A*belta+mu*omega)*(epsilon/(gamma+epsilon))*(2-l);
for i=1:length(l)
FVNj2(i)= (1+A*belta+mu*theta)/(1+A*belta+mu*omega)*(epsilon/(gamma+epsilon))*(2-l(i));
end
VNj2 = [FVNj2];


FGj2 =1/(1+sigma)*log(B)+theta/(1+sigma)*log(1/e*v*theta/(omega-theta)*(1+delta/R*(1+t)*l))+sigma/(1+sigma)*log(delta*(1+A*belta+mu*theta)/(1+A*belta+mu*omega)*(epsilon/(gamma+epsilon))*(2-l));
for i=1:length(l)
FGj2(i)=1/(1+sigma)*log(B)+theta/(1+sigma)*log(1/e*v*theta/(omega-theta)*(1+delta/R*(1+t)*l(i)))+sigma/(1+sigma)*log(delta*(1+A*belta+mu*theta)/(1+A*belta+mu*omega)*(epsilon/(gamma+epsilon))*(2-l(i)));
end

Hj2= [FGj2];

Gj2=exp(FGj2);

%sigma=0.1;gamma=0.1;
sigma=0.1;gamma=0.1;
l=0.4:0.01:0.6;%延迟退休到65；
FUj3=epsilon/(gamma+epsilon)*(1-l)-gamma/(gamma+epsilon);
for i=1:length(l)
FUj3(i)=epsilon/(gamma+epsilon)*(1-l(i))-gamma/(gamma+epsilon);
end
Uj3 = [FUj3];

FNj3=mu*(omega-theta)/(1+A*belta+mu*omega)*1/v*epsilon/(gamma+epsilon)*(2-l);
for i=1:length(l)
FNj3(i)=mu*(omega-theta)/(1+A*belta+mu*omega)*1/v*epsilon/(gamma+epsilon)*(2-l(i));
end
Nj3 = [FNj3];


FVNj3 = (1+A*belta+mu*theta)/(1+A*belta+mu*omega)*(epsilon/(gamma+epsilon))*(2-l);
for i=1:length(l)
FVNj3(i)= (1+A*belta+mu*theta)/(1+A*belta+mu*omega)*(epsilon/(gamma+epsilon))*(2-l(i));
end
VNj3 = [FVNj3];


FGj3 =1/(1+sigma)*log(B)+theta/(1+sigma)*log(1/e*v*theta/(omega-theta)*(1+delta/R*(1+t)*l))+sigma/(1+sigma)*log(delta*(1+A*belta+mu*theta)/(1+A*belta+mu*omega)*(epsilon/(gamma+epsilon))*(2-l));
for i=1:length(l)
FGj3(i)=1/(1+sigma)*log(B)+theta/(1+sigma)*log(1/e*v*theta/(omega-theta)*(1+delta/R*(1+t)*l(i)))+sigma/(1+sigma)*log(delta*(1+A*belta+mu*theta)/(1+A*belta+mu*omega)*(epsilon/(gamma+epsilon))*(2-l(i)));
end

Hj3= [FGj3];

Gj3=exp(FGj3);

%sigma=0.1;gamma=0.2;
sigma=0.1;gamma=0.2;
l=0.4:0.01:0.6;%延迟退休到65；
FUj5=epsilon/(gamma+epsilon)*(1-l)-gamma/(gamma+epsilon);
for i=1:length(l)
FUj5(i)=epsilon/(gamma+epsilon)*(1-l(i))-gamma/(gamma+epsilon);
end
Uj5 = [FUj5];

FNj5=mu*(omega-theta)/(1+A*belta+mu*omega)*1/v*epsilon/(gamma+epsilon)*(2-l);
for i=1:length(l)
FNj5(i)=mu*(omega-theta)/(1+A*belta+mu*omega)*1/v*epsilon/(gamma+epsilon)*(2-l(i));
end
Nj5 = [FNj5];


FVNj5 = (1+A*belta+mu*theta)/(1+A*belta+mu*omega)*(epsilon/(gamma+epsilon))*(2-l);
for i=1:length(l)
FVNj5(i)= (1+A*belta+mu*theta)/(1+A*belta+mu*omega)*(epsilon/(gamma+epsilon))*(2-l(i));
end
VNj5 = [FVNj5];


FGj5 =1/(1+sigma)*log(B)+theta/(1+sigma)*log(1/e*v*theta/(omega-theta)*(1+delta/R*(1+t)*l))+sigma/(1+sigma)*log(delta*(1+A*belta+mu*theta)/(1+A*belta+mu*omega)*(epsilon/(gamma+epsilon))*(2-l));
for i=1:length(l)
FGj5(i)=1/(1+sigma)*log(B)+theta/(1+sigma)*log(1/e*v*theta/(omega-theta)*(1+delta/R*(1+t)*l(i)))+sigma/(1+sigma)*log(delta*(1+A*belta+mu*theta)/(1+A*belta+mu*omega)*(epsilon/(gamma+epsilon))*(2-l(i)));
end

Hj5= [FGj5];

Gj5=exp(FGj5);

%%
%附图1
clear
theta=0.5;sigma=0.1;e=0.8;belta=0.33;gamma=0.15;omega=1;epsilon=0.4;alpha=0.25;delta=1;
n=1.14;g=2.5;mu=0.8;
l=0.4;
u=epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega);
syms R;
F=@(mu,theta,R,alpha,delta,l,e,belta,n,g,u)@(R)(1-R*mu*theta/(R+belta*R+(1-e)/e*mu*omega))*(1+u)-((1-R*mu*theta/(R+belta*R+(1-e)/e*mu*omega))+R*(omega-theta)*mu);
R=fzero(F(mu,theta,R,alpha,delta,l,e,belta,n,g,u),rand);

syms v
F2=@(R,omega,theta,mu,e,belta,n,v)@(v)R*(omega-theta)*mu/(((l-R*mu*theta/((belta*R*mu+(1-e)/e*mu*omega)+R))*(belta*R*mu+(1-e)/e*mu*omega+R)+R*(omega-theta)*mu))*(1+u)/v-n;
v=fzero(F2(R,omega,theta,mu,e,belta,n,v),rand);

B=log(g)*(1+sigma)-theta/(1+sigma)*log(1/e*theta*v/(omega-theta)*(1+l*delta/R))-sigma/(1+sigma)*log(delta*(1-v*n+u));

e=0.8;

l=1.38:0.01:1.8
plot(l,1/(1+sigma)*B+log(1/e*theta*v/(omega-theta)*(1+l*delta/R))*theta/(1+sigma)+log(delta*((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega))))*sigma/(1+sigma))
str='$$[\theta \delta(\varepsilon-\gamma)-\mu \varepsilon(1+r)] /[\varepsilon \delta(\mu+\theta)]$$'
text(1,1,str,'Interpreter','latex','FontSize',16)

set(gca,'linewidth',1.5); 
set(gca,'xtick',[]);
set(gca,'xticklabel',{});
set(gca,'YTick', []);
set(gca,'fontsize',12);
set (gca,'FontSize',12)
xlabel('退休年龄')
ylabel('青年人就业能力')


%%
%附图2
clear
theta=0.5;sigma=0.1;e=0.8;belta=0.33;gamma=0.15;omega=1;epsilon=0.4;alpha=0.25;delta=1;
n=1.14;g=2.5;mu=0.8;
l=0.4;
u=epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega);
syms R;
F=@(mu,theta,R,alpha,delta,l,e,belta,n,g,u)@(R)(1-R*mu*theta/(R+belta*R+(1-e)/e*mu*omega))*(1+u)-((1-R*mu*theta/(R+belta*R+(1-e)/e*mu*omega))+R*(omega-theta)*mu);
R=fzero(F(mu,theta,R,alpha,delta,l,e,belta,n,g,u),rand);

syms v
F2=@(R,omega,theta,mu,e,belta,n,v)@(v)R*(omega-theta)*mu/(((l-R*mu*theta/((belta*R*mu+(1-e)/e*mu*omega)+R))*(belta*R*mu+(1-e)/e*mu*omega+R)+R*(omega-theta)*mu))*(1+u)/v-n;
v=fzero(F2(R,omega,theta,mu,e,belta,n,v),rand);

B=log(g)*(1+sigma)-theta/(1+sigma)*log(1/e*theta*v/(omega-theta)*(1+l*delta/R))-sigma/(1+sigma)*log(delta*(1-v*n+u));

e=0.8;

l=0.4:0.1:1.8
FG=1/(1+sigma)*B+log(1/e*theta*v/(omega-theta)*(1+l*delta/R))*theta/(1+sigma)+log(delta*((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega))))*sigma/(1+sigma);
for i=1:length(l)
FG(i)=1/(1+sigma)*B+log(1/e*theta*v/(omega-theta)*(1+l(i)*delta/R))*theta/(1+sigma)+log(delta*((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega))))*sigma/(1+sigma);
end
G = exp(FG);

l=0.4:0.01:0.6
FG1=1/(1+sigma)*B+log(1/e*theta*v/(omega-theta)*(1+l*delta/R))*theta/(1+sigma)+log(delta*((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l)-gamma/(gamma+epsilon*omega))))*sigma/(1+sigma);
for i=1:length(l)
FG1(i)=1/(1+sigma)*B+log(1/e*theta*v/(omega-theta)*(1+l(i)*delta/R))*theta/(1+sigma)+log(delta*((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)/((1-R*mu*theta/(belta*R+(1-e)/e*mu*omega+R))*(belta*R+(1-e)/e*mu*omega+R)+(omega-theta)*R*mu)*(1+epsilon*omega/(epsilon*omega+gamma)*(1-l(i))-gamma/(gamma+epsilon*omega))))*sigma/(1+sigma);
end

G1 = exp(FG1);

clf;
l1=0.4:0.01:0.6;
subplot(1,2,1)
plot(l1,G1);
l2=0.4:0.1:1.8;
subplot(1,2,2)
plot(l2,G);

subplot(1,2,1);
set(gca,'linewidth',1.5); 
set(gca,'xtick',[]);
set(gca,'xticklabel',{});
set(gca,'fontsize',12);
set (gca,'FontSize',12)
xlabel('退休年龄')
ylabel('青年人就业能力')

subplot(1,2,2);
set(gca,'linewidth',1.5); 
set(gca,'xtick',[]);
set(gca,'xticklabel',{});
set(gca,'fontsize',12);
set (gca,'FontSize',12)
xlabel('退休年龄')
ylabel('青年人就业能力')

