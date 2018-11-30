function f=plotGravity(G1,G2,G3,G4)
factor=6;
int=1/factor;
[a,b]=size(G1);
tend=a/factor;
tmax=100;
g1=nanmean(G1,2);
g2=nanmean(G2,2);
g3=nanmean(G3,2);
g4=nanmean(G4,2);
maxy=max([g1;g2;g3;g4]);
tvec=(int:int:tend);
fs=15; lw=2;
figure
plot(tvec,g1,'-','linewidth',lw)
hold on
plot(tvec,g2,'-','linewidth',lw)
plot(tvec,g3,'-','linewidth',lw)
plot(tvec,g4,'-','linewidth',lw)
legend('PG','OG','PR','OR')
xlabel('Time (days)','FontSize',fs)
ylabel('Prevalence','FontSize',fs)
axis ([0,tmax,0,maxy]);
%axis tight
set(gca,'FontSize',fs);
grid on
grid minor
box on
