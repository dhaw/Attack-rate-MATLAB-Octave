function f=plotAR(F,NN,a,b)%,tit)
%Trim:
%{
NN=trimbyk(NN,a,b); n=length(NN); logN=log10(NN);
num=size(F,2); F2=zeros(n,num);
%{
for i=1:num
    F2(:,i)=trimbyk(F(:,i),a,b);
end
F=F2;
%}
meanA=nansum(F2(:,1).*NN)./sum(NN);
%}
n=length(NN); logN=log10(NN);
%For CAR/PI:
A=F;
AA=[min(A,[],2),max(A,[],2)];
AA=prctile(A,[25,75],2);
A=nanmean(A,2);
howmany=min(n,100);
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
meanA=sum(A.*NN)./sum(Anan.*NN);

Ac=A; Ac(logN==-inf)=[]; logNc=logN; logNc(logN==-inf)=[];
cc=corrcoef(logNc,Ac); cc=cc(2);
f=[meanA,cc];

figure
%suptitle(tit)
%fs=40; ms=30; lwx=3; lw=3; %For presentations
fs=12; ms=5;%12; 20 %fs=18 for paper figures
lwx=1; lw=1.5; %lwx=1.5
col1=[.165,.31,.431];%[1,0,0];%[0.0512,0.4600,0.8633];%[.165,.31,.431];
col2=[.447,.553,.647];
col3=[0,0,0];
%hold on
semilogx(NN,A,'o','color',col2,'markersize',ms,'LineWidth',lwx,'markerfacecolor',col2);
%plot(NN,A,'o','color',col1,'markersize',ms,'LineWidth',lwx);
hold on
semilogx([1,maxN],[meanA,meanA],'k--','linewidth',lw);

%errorbar(Nsam,meansam,AA(:,1),AA(:,2),'ko','markerfacecolor','w','markeredgecolor','k','markersize',ms);

%semilogx([1,maxN],[0,0],'k-','linewidth',lw);
%plot([0,maxN],[meanA,meanA],'k--','linewidth',lw);
%CAR/PI:
Aall=F;%[Arur;Aurb];
maxY=max(max(Aall)); minY=min(min(Aall));
%Comment out to remove error bars:
%{
for i=1:length(Nsam)
    semilogx([Nsam(i),Nsam(i)],Asam(i,:),'color','k','linewidth',lw);%Plotted over above [.071,.208,.322]
    semilogx(Nsam(i),meansam(i),'o','color','k','markerfacecolor','w','markeredgecolor','k','markersize',ms);
end
%}
xlabel('Pop. density','FontSize',fs);
ylabel('Attack rate')
%xlabel('Proportion aged 65+','FontSize',fs);
%ylabel(strcat('Attack rate (',tit,')'),'FontSize',fs);
set(gca,'FontSize',fs,'xtick',[1,10,100,1000,10000,100000]);
minY=max(minY-.1,0); maxY=min(maxY+.1,1);
axis([0,maxN,0,maxY])%minY,maxY])
%set(gca,'yticklabels','')
grid on
grid minor
box on
hold off
end