function [f,g]=sensAnalAR2(C1,C2,C3,C4)%(C1,C2,C3,C4)(C,Q)
aa=(.1:.1:4);
laa=length(aa);
X=zeros(5,laa);
Y=zeros(4,laa);
celldist=1;
%%Within-age-group local attack rates - diff from global means
%
for i=1:laa
    aai=aa(i);
    [gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,~]=prepFluAge(C1,C2,C3,C4,1.8,0,0,celldist,aai); beta3=1;
    %[gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,~]=prepFluAgeLocs(C,Q,0,1,aai);
    [f,~,~]=finalSizeMultiAgeOut(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,1,2,.00001,1);
    F=f.*NNrep;%Local absolute attack
    Fs=reshape(F,n,4); Fs=sum(Fs,1); Fs=Fs';
    Ns=reshape(NNbar,n,4); Ns=sum(Ns,1); Ns=Ns'; 
    globmeans=Fs./Ns;
    F=F./NNbar;%Within-age-group local attack rates
    F(NNbar==0)=0;
    F=F-kron(globmeans,ones(n,1));
    X(:,i)=prctile(F,[0,25,50,75,100]);
    Y(:,i)=globmeans;
end
figure
f=X;
g=Y;
minx=min(min(X)); maxx=max(max(X));
plot([.5,laa+.5],[0,0],'k-','linewidth',1)
hold on
h=boxplot(X,'whisker',10);
set(h,{'linew'},{1})
set(gca,'xtick',(5:5:laa),'xticklabel',aa(5:5:end));
hold off
set(gca,'fontsize',12)
%set(gca,'xtick',0:6)
axis ([.5,laa+.5,minx-.05,maxx+.05])
%title('Before')
xlabel('Offset a')
ylabel('Attack rate (diff)')
grid on
grid minor
box on
%}
%%Global attack
%{
for i=1:laa
    aai=aa(i);
    [gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,~]=prepFluAge(C1,C2,C3,C4,1.8,0,0,celldist,aai); beta3=1;
    %[gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,~]=prepFluAgeLocs(C,Q,0,1,aai);
    [f,~,~]=finalSizeMultiAgeOut(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,1,2,.00001,1);
    F=f.*NNrep;%Local absolute attack
    F=reshape(F,n,4);
    F=sum(F,2)./NN;
    F(NN==0)=0;
    X(:,i)=prctile(F,[0,25,50,75,100]);
end
f=X;
g=Y;
minx=min(min(X)); maxx=max(max(X));
figure
hold on
h=boxplot(X,'whisker',10);
set(h,{'linew'},{1})
set(gca,'xtick',(5:5:laa),'xticklabel',aa(5:5:end));
hold off
set(gca,'fontsize',12)
%set(gca,'xtick',0:6)
axis ([.5,laa+.5,0,maxx+.05])
%title('Before')
xlabel('Offset a')
ylabel('Attack rate')
grid on
grid minor
box on
%}