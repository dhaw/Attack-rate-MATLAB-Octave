function f=plotARfscape(F,NN,Is)%,tit)
meanA=nansum(F(:,1).*NN)./sum(NN);
%}
n=length(NN); logN=log10(NN);
%For CAR/PI:
A=F;
maxN=max(NN);%logN);
Anan=1-isnan(A);
A(Anan==0)=0;
meanA=sum(A.*NN)./sum(Anan.*NN);

Ac=A; Ac(logN==-inf)=[]; logNc=logN; logNc(logN==-inf)=[];
cc=corrcoef(logNc,Ac); cc=cc(2);
f=[meanA,cc];

fscapeinds=find(Is);
locs=Is(fscapeinds);
NB=NN(locs);
B=A(locs,:);


figure
%suptitle(tit)
%fs=40; ms=30; lwx=3; lw=3; %For presentations
fs=18; ms=5;%12; 20
lwx=1; lw=1.5; %lwx=1.5
col1=[.165,.31,.431];%[1,0,0];%[0.0512,0.4600,0.8633];%[.165,.31,.431];
col2=[.447,.553,.647];
col3=[0,0,0];
%hold on
semilogx(NN,A,'o','color',col2,'markersize',ms,'LineWidth',lwx);
%plot(NN,A,'o','color',col1,'markersize',ms,'LineWidth',lwx);
hold on
semilogx(NB,B,'o','color',[1,0,0],'markerfacecolor',[1,0,0],'markersize',ms,'LineWidth',lwx);
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
set(gca,'FontSize',fs);
minY=max(minY-.1,0); maxY=min(maxY+.1,1);
axis([0,maxN,0,maxY])%minY,maxY])
%set(gca,'yticklabels','')
grid on
grid minor
box on
hold off
end