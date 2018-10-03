function f=plotMAnospace(Asus,Arec,NN,ar)
fs=15; lw=2;%30;%Font size
tauend=size(Asus,2);
%
figure
T=1:tauend;
hold on
if ar==1
    plot(T,Arec,'-','linewidth',lw)%,'color','k')%'o-'
else
    plot(T,Asus,'-','linewidth',lw)%,'color','k')%'o:'
end

xlabel('Time (years)','FontSize',fs)
%ylabel('Proportion immune','FontSize',fs)
if ar==1
    ylabel('Prop. newly recovered','FontSize',fs)
else
    ylabel('Prop. susceptible','FontSize',fs)
end
set(gca,'FontSize',fs);
%maxY=max(max([Acum;A])); maxY=min(maxY+.1,1);
axis tight%([0,tauend,0,1]);%maxY])

if ar==1
    legend('R^1','R^2')
else
    legend('S_{11}','S_{10}','S_{01}','S_{00}')
end
grid on
grid minor
hold off