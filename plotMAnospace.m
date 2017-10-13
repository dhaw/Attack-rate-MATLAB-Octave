function f=plotMAnospace(A,Acum,NN)
fs=15;%30;%Font size
tauend=size(A,2);
%
figure
T=1:tauend;
hold on
plot(T,Acum,'o-','linewidth',1.5,'color','k')
plot(T,A,'o:','linewidth',1.5,'color','k')
xlabel('Time (years)','FontSize',fs)
ylabel('Proportion immune','FontSize',fs)
set(gca,'FontSize',fs);
maxY=max(max([Acum;A])); maxY=min(maxY+.1,1);
axis([1,tauend,0,maxY])
grid on
hold off