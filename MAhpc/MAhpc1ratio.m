function f=MAhpc1ratio
load('forMAhpc.mat')
eps=(0:.1:1);
leps=length(eps);
%chi=(0:.01:1);
thresh=10^(-4);
%
thismany=5;%Random ICs in fSMR - loop length
tauend=200;
burn=100;
lt=tauend-burn;
%years=tauend-burn;
isdual=3;
solvetype=2;
numseed=10^(-8);
[gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,ages0]=prepFluAgeLocs(C,Qeven,0,1);
%X=zeros(n,lt,thismany);
c=cell(1,thismany);
ceps=cell(1,leps);
miny=10; maxy=0;
for j=1:leps%thismany
    epsj=eps(j);
    %try
    [f,g]=finalSizeMulti(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,isdual,solvetype,numseed,epsj);
    x=zeros(n,lt);
    for i=1:tauend-1
        x(:,i)=f(:,i+1)./f(:,i); 
    end
    %c{j}=unique(x);
    cepsj=unique(x);
    ceps{j}=cepsj;
    miny=min(miny,min(cepsj));
    maxy=max(maxy,max(cepsj));
end
f=ceps;

fs=12; lw=2; ms=7.5;
figure
hold on
for j=1:leps;
    plot(eps(j),ceps{j},'.','markersize',ms,'color',[0,0,0])
end
hold off
%axis([.5,thismany+.5,miny,maxy])
axis([0,1,miny,maxy])
xlabel('Trial')
ylabel('Ratio')
set(gca,'fontsize',fs)
%legend('Total immune','Newly infected')
grid on
grid minor
box on