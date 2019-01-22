function [f,g]=ar2ic(C,Q)%(X)(C1,C2,C3,C4)(C,Q)
eps=(0:.1:1);
leps=length(eps);
%
X=zeros(5,leps);
Y=zeros(1,leps);
beta3=1;
%[gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,ages0]=prepFluAge(C1,C2,C3,C4,1.8,0,1);
[gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,ages0]=prepFluAgeLocs(C,Q,0,1);
NN4=NNbar(3*n+1:end); NN0=NN; NN0(NN==0)=1; icmax4=NN4./NN0;
NN3=NNbar(2*n+1:3*n); icmax3=NN3./NN0;
for i=1:leps
    epsi=eps(i);
    [f,~,~]=finalSizeMulti(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,1,2,.00001,epsi);
    ic=epsi*icmax4;
    %ic=epsi*(icmax4+12/45*icmax3);
    f=f-ic;%Attack only
    xi=prctile(f,[0,25,50,75,100]);
    X(:,i)=xi;
    glob=sum(NN.*f)/sum(NN);
    Y(i)=glob;
end
f=X;
g=Y;
%}
figure;
%plot(0:6,HH(:,1),'k-','linewidth',1)
plot((1:2:leps),Y(1:2:leps),'k-','linewidth',1.5)
hold on
h=boxplot(X,'whisker',10);%,'positions',res,'labels',res);
set(h,{'linew'},{1})
set(gca,'xtick',(1:2:leps),'xticklabel',eps(1:2:leps));
hold off
set(gca,'fontsize',12)
%set(gca,'xtick',0:6)
minx=min(min(X)); maxx=max(max(X));
axis ([.5,leps+.5,0,1])%,maxx+.05])%max(minx-.05,0)
%title('Before')
xlabel('\epsilon')
ylabel('Attack rate')
grid on
grid minor
box on