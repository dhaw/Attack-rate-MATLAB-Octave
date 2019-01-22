function f=plotARsplit(F,NN,NNprop,prop,a,b,tit)
%NN - prop of age group - x axis
%NNprop - population to split
%prop - proportion in population

%Single year only
%Trim:
%
NN=trimbyk(NN,a,b); n=length(NN); logN=log10(NN);
NNprop=trimbyk(NNprop,a,b);
num=size(F,2); F2=zeros(n,num);
A=trimbyk(F,a,b);
%}
n=length(NN); logN=log10(NN);
%
maxN=max(NN);%logN);
Anan=1-isnan(A);
A(Anan==0)=0;
meanA=sum(A.*NN)./sum(Anan.*NN);
Ac=A; Ac(logN==-inf)=[]; logNc=logN; logNc(logN==-inf)=[];
cc=corrcoef(logNc,Ac); cc=cc(2);
f=[meanA,cc];

propG=A; propG(NNprop<prop)=-1;
propL=A; propL(NNprop>prop)=-1;

figure
%suptitle(tit)
%fs=40; ms=30; lwx=3; lw=3; %For presentations
fs=12; ms=5;%12; 20
lwx=1; lw=1.5; %lwx=1.5
col1=[0,0,.5];%[0,0,1];%[.8,0,0];%[0.0512,0.4600,0.8633];%[.165,.31,.431];
col2=[.7,.7,.7];%[.3,.3,.3];%[0.4612,0.3628,0.1509];%[.447,.553,.647];
%hold on
%semilogx(NN,A,'o','color',col1,'markersize',ms,'LineWidth',lwx);
plot(NN,propL,'o','color',col2,'markersize',ms,'LineWidth',lwx,'markerfacecolor',col2);
hold on
plot(NN,propG,'o','color',col1,'markersize',ms,'LineWidth',lwx);

%semilogx([1,maxN],[meanA,meanA],'k--','linewidth',lw);
%semilogx([1,maxN],[0,0],'k-','linewidth',lw);
plot([0,maxN],[meanA,meanA],'--','linewidth',lw,'color',[.5,0,0]);
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
%xlabel('Population density (N/cell)','FontSize',fs);
xlabel('Proportion aged 0-4','FontSize',fs);
%xlabel('Proportion aged 5-19','FontSize',fs);
%xlabel('Proportion aged 20-64','FontSize',fs);
%xlabel('Proportion aged 65+','FontSize',fs);
%
%ylabel(strcat('Attack rate (',tit,')'),'FontSize',fs);
ylabel('Attack rate')
set(gca,'FontSize',fs);
minY=max(minY-.1,0); maxY=min(maxY+.1,1);
axis ([0,maxN,0,maxY])%minY,maxY])
%legend('0-4<av.','0-4>av.','location','NW')
%legend('5-19<av.','5-19>av.','location','NW')
legend('20-64<av.','20-64>av.','location','NW')
%legend('65+<av.','65+>av.','location','NW')
grid on
grid minor
hold off
end