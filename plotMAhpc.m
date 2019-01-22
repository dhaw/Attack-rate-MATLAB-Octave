function f=plotMAhpc(X,Y,Z)
tauend=size(Z,1);
timeint=ceil(tauend/5);
thismany=size(Y,2);
trialint=ceil(thismany/5);

fs=15; lw=2;
%%Variance:
%{
figure
boxplot(Y)
xlabel('Trial','FontSize',fs)
ylabel('Variance of rank','FontSize',fs)
%axis ([.5,thismany+.5,-1,1]);
axis tight
set(gca,'FontSize',fs,'xtick',0:trialint:tauend);
grid on
grid minor
box on
%}
%%Z in time:
%{
figure
plot(1:tauend,Z,'-','linewidth',lw)%,'color',[.5,0,0])%'o-'
xlabel('Time (years)','FontSize',fs)
ylabel('Correlation with v_1','FontSize',fs)
axis tight%([0,tauend,-1,1]);
%axis tight
set(gca,'FontSize',fs,'xtick',0:timeint:tauend);
grid on
grid minor
box on
%}
%%Histogram of cc's:
[a,b]=size(Z);
Z2=reshape(Z,a*b,1);
%Z2(Z==0)=[];
fs=12;
lw=2;
figure
hist(Z2,20,'facecolor',[0,0,.5]);
xlabel('Corr. coef.')
ylabel('Frequency')
set(gca,'fontsize',fs)
axis([-1,1,0,max(hist(Z2))])
grid on
grid minor
box on

%%
%{
figure
plot(1:tauend-1,x,'linewidth',1)
xlabel('Time (years)')
ylabel('Ratio')
set(gca,'fontsize',15)
grid on
grid minor
box on;
%}