function [f,g,h]=MAhpcGeom
eps=(0:.1:1);
leps=length(eps);
thismany=2;%Random ICs in fSMR - loop length
isdual=0;
solvetype=2;
numseed=10^(-8);
%%Fluscape (cells only):
%{
load('fluscape.mat')
[Df,xf,yf,Ds,Is]=fluscapeOnlyLocalAverage(lscan,fluscapeLocations);
fscapeind=find(Is);
[lscanNew,r]=fluscapeNNr(Df,Ds);
[gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,ages0]=prepFluAgeLocsFscape(lscanNew(fscapeind),r(fscapeind,fscapeind),0,1);
%}
%%Pittsburgh
load('forMAhpc.mat')
[gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,ages0]=prepFluAgeLocs(C,Qeven,0,1);
%Pandemic case:
randic=0;
[fp,gp,Dp]=finalSizeMulti(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,isdual,solvetype,numseed,1,randic,1);
xi=[(1:n)',fp];
xi=sortrows(xi,2);
prank=xi(:,1);
%
X=zeros(n,leps,thismany);
Y=zeros(leps,thismany);
Z=Y;
randic=1;
tauend=100;
for i=1:leps
    epsi=eps(i);
    if epsi==0
        tauend=500;
    end
    for j=1:thismany
        [f,g,D]=finalSizeMulti(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,isdual,solvetype,numseed,epsi,randic,tauend);
        %
        fx=f(:,end);
        X(:,i,j)=fx;
        ccij=corrcoef(fx,fp);
        Y(i,j)=ccij(1,2);
        xi=[(1:n)',fx];
        xi=sortrows(xi,2);
        xi=xi(:,1);
        ccij=corrcoef(xi,prank);
        Z(i,j)=ccij(1,2);
    end
end
f=X;
g=Y;
h=Z;
%save('MA1geom','X','Y')
end