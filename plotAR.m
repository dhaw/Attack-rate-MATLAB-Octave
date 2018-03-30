<<<<<<< HEAD
function f=plotAR(F,NN,a,b,tit)
%Trim:
%
NN=trimbyk(NN,a,b); n=length(NN); logN=log10(NN);
num=size(F,2); F2=zeros(n,num);
for i=1:num
    F2(:,i)=trimbyk(F(:,i),a,b);
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
meanA=sum(A.*NN)./sum(Anan.*NN);

Ac=A; Ac(logN==-inf)=[]; logNc=logN; logNc(logN==-inf)=[];
cc=corrcoef(logNc,Ac); cc=cc(2);
f=[meanA,cc];

figure
%suptitle(tit)
%fs=40; ms=30; lwx=3; lw=3; %For presentations
fs=12; ms=5;%12; 20
lwx=1; lw=1.5; %lwx=1.5
col1=[1,0,0];%[0.0512,0.4600,0.8633];%[165,.31,.431];
col2=[.447,.553,.647];
%hold on
semilogx(NN,A,'o','color',col1,'markersize',ms,'LineWidth',lwx);
%plot(NN,A,'o','color',col1,'markersize',ms,'LineWidth',lwx);
hold on
semilogx([1,maxN],[meanA,meanA],'k--','linewidth',lw);
%semilogx([1,maxN],[0,0],'k-','linewidth',lw);
%plot([0,maxN],[meanA,meanA],'k--','linewidth',lw);
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
xlabel('Population density (N/cell)','FontSize',fs);
%xlabel('Proportion aged 65+','FontSize',fs);
ylabel(strcat('Attack rate (',tit,')'),'FontSize',fs);
set(gca,'FontSize',fs);
minY=max(minY-.1,0); maxY=min(maxY+.1,1);
axis ([0,maxN,0,maxY])%minY,maxY])
grid on
grid minor
hold off
=======
function f=plotAR(F,NN,a,b,tit)
%Trim:
%
NN=trimbyk(NN,a,b); n=length(NN); logN=log10(NN);
num=size(F,2); F2=zeros(n,num);
for i=1:num
    F2(:,i)=trimbyk(F(:,i),a,b);
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
meanA=sum(A.*NN)./sum(Anan.*NN);

Ac=A; Ac(logN==-inf)=[]; logNc=logN; logNc(logN==-inf)=[];
cc=corrcoef(logNc,Ac); cc=cc(2);
f=[meanA,cc];

figure
%suptitle(tit)
%fs=40; ms=30; lwx=3; lw=3; %For presentations
fs=12; ms=5;%12; 20
lwx=1; lw=1.5; %lwx=1.5
col1=[0.0512,0.4600,0.8633];%[165,.31,.431];
col2=[.447,.553,.647];
%hold on
semilogx(NN,A,'o','color',col1,'markersize',ms,'LineWidth',lwx);
%plot(NN,A,'o','color',col1,'markersize',ms,'LineWidth',lwx);
hold on
semilogx([1,maxN],[meanA,meanA],'k--','linewidth',lw);
%semilogx([1,maxN],[0,0],'k-','linewidth',lw);
%plot([0,maxN],[meanA,meanA],'k--','linewidth',lw);
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
xlabel('Population density (N/cell)','FontSize',fs);
%xlabel('Proportion aged 65+','FontSize',fs);
ylabel(strcat('Attack rate (',tit,')'),'FontSize',fs);
set(gca,'FontSize',fs);
minY=max(minY-.1,0); maxY=min(maxY+.1,1);
axis ([0,maxN,0,maxY])%minY,maxY])
grid on
grid minor
hold off
>>>>>>> acfad3e4da850e525286e08138beef9f06a8c2cc
end