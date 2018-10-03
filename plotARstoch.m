function f=plotARstoch(F,NN,a,b)%,tit)
%Trim:
%
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
%AA=[min(A,[],2),max(A,[],2)];
AA=prctile(A,[25,50,75],2);
A=nanmean(A,2);
AA=abs(repmat(AA(:,2),1,2)-AA(:,[1,3]));

howmany=100;
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

figure
fs=18; ms=5;
lwx=1; lw=1.5;
col1=[.165,.31,.431];
col2=[.447,.553,.647];
col3=[0,0,0];
%hold on
semilogx(NN,A,'o','color',col2,'markersize',ms,'LineWidth',lwx);
%plot(NN,A,'o','color',col1,'markersize',ms,'LineWidth',lwx);
hold on
%semilogx([1,maxN],[meanA,meanA],'k--','linewidth',lw);
errorbar(Nsam,meansam,Asam(:,1),Asam(:,2),'ko','markerfacecolor','w','markeredgecolor','k','markersize',ms,'LineWidth',lwx);

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
set(gca,'FontSize',fs,'xtick',[1,10,100,1000,10000,100000]);
set(gca,'FontSize',fs);
minY=max(minY-.1,0); maxY=min(maxY+.1,1);
axis tight%([0,maxN,0,maxY])%minY,maxY])
grid on
grid minor
box on
hold off
end