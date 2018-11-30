function f=plotARscatter(NN,f0,f1,f2)%S/D/I - as "isdual"
fs=12; ms=5;%12; 20
lwx=1; lw=1.5; %lwx=1.5
col1=[.165,.31,.431];%[1,0,0];%[0.0512,0.4600,0.8633];%[.165,.31,.431];
col2=[.447,.553,.647];
col3=[0,0,0];
figure
maxN=max(NN);
semilogx([1,maxN],[0,0],'k-','linewidth',lw);
hold on
%h1=semilogx(NN,f1,'o','color',col3,'markerfacecolor',col3,'markersize',ms,'LineWidth',lwx);
h2=semilogx(NN,f0,'o','color',col1,'markerfacecolor',col1,'markersize',ms,'LineWidth',lwx);
%h3=semilogx(NN,f2,'o','color',col2,'markerfacecolor',col2,'markersize',ms,'LineWidth',lwx);

%h2=semilogx(NN,f1-f0,'o','color',col1,'markerfacecolor',col1,'markersize',ms,'LineWidth',lwx);
%h3=semilogx(NN,f1-f2,'o','color',col2,'markerfacecolor',col2,'markersize',ms,'LineWidth',lwx);

%=[min(Rs0),max([R00;R03;R06])];
%plot(v,v,'k-','linewidth',2)
xlabel('Pop. density','FontSize',fs);
%xlabel('R_0','FontSize',fs);
ylabel('Attack rate','FontSize',fs);
set(gca,'FontSize',fs);
axis tight

%legend([h1,h2,h3],'DM','SM','IM','location','NW')
%legend([h2,h3],'SM','IM','location','NW')
grid on
grid minor
box on
hold off