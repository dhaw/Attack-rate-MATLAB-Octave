function f=plotMAhpc2sub(X,Y,Z)%H1,H3,prop H1

[a,b,c]=size(X);

%Concatenate simulations:
%X=reshape(X,a,b*c,1);

%Distribution of variances:
X=nanvar(X,[],2);
X=reshape(X,[a,c,1]);
Y=nanvar(Y,[],2);
Y=reshape(Y,[a,c,1]);
Z=nanvar(Z,[],2);
Z=reshape(Z,[a,c,1]);
%X=sqrt(X);
%Y=sqrt(Y);
%Z=sqrt(Z);
eps=(0:.05:1);
fs=12; lw=2;
cmod=.3;
figure;
hold on
plot(eps,X,'-','linewidth',lw)
plot(eps,Y,'-','linewidth',lw)
plot(eps,Z,'--','color',cmod*[1,1,1],'linewidth',lw)
hold off
hold off
set(gca,'fontsize',fs)
axis ([0,eps(end),0,.2])
xlabel('\epsilon')
ylabel('\sigma^2(cc_{pand}), \chi=0.8')
%ylabel('\sigma(cc_{pand}), \chi=0.2')
legend('H1N1','H3N2','Total','location','NE')
grid on
grid minor
box on
%{
leps=length(eps);
xbox=prctile(X,[0,25,50,75,100],2);
ybox=prctile(Y,[0,25,50,75,100],2);
zbox=prctile(Z,[0,25,50,75,100],2);
subplot(1,3,1)
%hold on
%plot(1:leps,y,'k-','linewidth',lw)
%plot(1:leps,z,'k--','linewidth',lw)%,'color',[0,.5,0]
h=boxplot(xbox','whisker',10);
set(h,{'linew'},{1})
set(gca,'xtick',(1:2:leps),'xticklabel',eps(1:2:leps));
hold off
set(gca,'fontsize',fs)
axis ([0,leps+1,-1,1])
xlabel('\epsilon')
%ylabel('\sigma(cc_{pand}), \chi=0.4')
ylabel('\sigma(cc_{pand}), \chi=0.2')
grid on
grid minor
box on
%
subplot(1,3,2)
%hold on
%plot(1:leps,y,'k-','linewidth',lw)
%plot(1:leps,z,'k--','linewidth',lw)%,'color',[0,.5,0]
h=boxplot(ybox','whisker',10);
set(h,{'linew'},{1})
set(gca,'xtick',(1:2:leps),'xticklabel',eps(1:2:leps));
hold off
set(gca,'fontsize',fs)
axis ([0,leps+1,-1,1])
xlabel('\epsilon')
%ylabel('Correlation')
grid on
grid minor
box on
%
subplot(1,3,3)
%hold on
%plot(1:leps,y,'k-','linewidth',lw)
%plot(1:leps,z,'k--','linewidth',lw)%,'color',[0,.5,0]
h=boxplot(zbox','whisker',10);
set(h,{'linew'},{1})
set(gca,'xtick',(1:2:leps),'xticklabel',eps(1:2:leps));
hold off
set(gca,'fontsize',fs)
axis ([0,leps+1,-1,1])
xlabel('\epsilon')
%ylabel('Correlation')
grid on
grid minor
box on
%}
%%Drift - rank order and epsilon:
%

%}