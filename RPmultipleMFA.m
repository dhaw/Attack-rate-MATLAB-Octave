function f=RPmultipleMFA
%Change epiRP to output g=Ysum
p=(0:.5:6);
lp=length(p);
tend=360;
plotto=180;
kangle=10;
h=4;
maxY=0;

figure
fs=12; lw=2;
cmap=parula(lp);

[gamma,NN,n,K,L,Vi,Viover,kangle,h,w,rhohat,isflat,beta,boxLat,boxLong]=RPprep(1,3,kangle,h,p(1));
loc=randsample(n,10);
[tout,Ysum]=RPepi(gamma,NN,n,K,L,Vi,Viover,beta,kangle,h,w,rhohat,isflat,10^(-4),1,boxLat,boxLong,loc);%Sample=1 - redundant
semilogy(tout,Ysum/n/h,'color',cmap(1,:),'linewidth',lw);
maxY=max(maxY,max(Ysum));

hold on
for i=2:lp
    [gamma,NN,n,K,L,Vi,Viover,kangle,h,w,rhohat,isflat,beta,boxLat,boxLong]=RPprep(1,3,kangle,h,p(i));
    [tout,Ysum]=RPepi(gamma,NN,n,K,L,Vi,Viover,beta,kangle,h,w,rhohat,isflat,10^(-4),1,boxLat,boxLong,loc);
    semilogy(tout,Ysum/n/h,'color',cmap(i,:),'linewidth',lw);
    maxY=max(maxY,max(Ysum));
end
xlabel('Time (days)','FontSize',fs);
ylabel('Incidence','FontSize',fs);%Prevalence
set(gca,'FontSize',fs);
axis ([30,min(tend,plotto),10^(-4),maxY/n/h])
legend('\alpha=2','\alpha=2.5','\alpha=3','\alpha=3.5','\alpha=4','\alpha=4.5','\alpha=5','\alpha=5.5','\alpha=6','location','SE')
grid on
grid minor
hold off