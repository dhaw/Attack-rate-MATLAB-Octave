function f=plotAR(F,NN)
%Trim:
%{
NN=trimbyk(NN); n=length(NN); logN=log10(NN);
num=size(F,2); F2=zeros(n,num);
for i=1:num
    F2(:,i)=trimbyk(F(:,i));
end
F=F2;
meanA=nansum(F2(:,1).*NN)./sum(NN);
%}
n=length(NN); logN=log10(NN);
%For CAR/PI:
A=F;
AA=[min(A,[],2),max(A,[],2)];
A=nanmean(A,2);
howmany=1;
isam=randsample(n,howmany);
Asam=zeros(howmany,2);
Nsam=zeros(howmany,1);
meansam=Nsam;
for i=1:howmany
        Asam(i,:)=AA(isam(i),:);
        Nsam(i,1)=NN(isam(i));
        meansam(i,1)=median(F(isam(i),:),2);
end
%
maxN=max(NN);%logN);
Anan=1-isnan(A);
A(Anan==0)=0;
%meanA=sum(A.*NN)./sum(Anan.*NN);

Ac=A; Ac(logN==-inf)=[]; logNc=logN; logNc(logN==-inf)=[];
cc=corrcoef(logNc,Ac); cc=cc(2);
f=[meanA,cc];

figure
fs=40;%Font size
%hold on
semilogx(NN,A,'x','color',[.165,.31,.431],'markersize',30,'LineWidth',3);%[.165,.31,.431][.447,.553,.647]
hold on
semilogx([1,maxN],[meanA,meanA],'k--','linewidth',3);
%CAR/PI:
Aall=F;%[Arur;Aurb];
maxY=max(max(Aall)); minY=min(min(Aall));
%Comment out to remove error bars:
%{
for i=1:length(Nsam)
    semilogx([Nsam(i),Nsam(i)],Asam(i,:),'color',[.071,.208,.322],'linewidth',3);%Plotted over above
    semilogx(Nsam(i),meansam(i),'o','color',[.071,.208,.322],'markerfacecolor',[.071,.208,.322],'markeredgecolor',[.071,.208,.322],'markersize',12);
end
%}
xlabel('Population density (N/km^2)','FontSize',fs);
ylabel('Proportion immune','FontSize',fs);
set(gca,'FontSize',fs);
minY=max(minY-.1,0); maxY=min(maxY+.1,1);
axis([0,maxN,0,maxY])%minY,maxY])
grid on
hold off
end