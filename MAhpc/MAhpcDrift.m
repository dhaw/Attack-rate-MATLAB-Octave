function [f,g]=MAhpcDrift%(lscan,fluscapeLocations)
thismany=2;%Random ICs in fSMR - loop length
tauend=200;
burn=100;
years=tauend-burn;
isdual=0;
solvetype=2;
numseed=10^(-5);
eps=(0:.1:1);
leps=length(eps);
%
load('forMAhpc.mat')
[gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,ages0]=prepFluAgeLocs(C,Qeven,0,1);
%}
%{
load('fluscape.mat')
[Df,xf,yf,Ds,Is]=fluscapeOnlyLocalAverage(lscan,fluscapeLocations);
fscapeind=find(Is);
[lscanNew,r]=fluscapeNNr(Df,Ds);
[gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,ages0]=prepFluAgeLocsFscape(lscanNew(fscapeind),r(fscapeind,fscapeind),0,1);
%}
[fp,gp,~]=finalSizeMulti(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,isdual,solvetype,numseed,0,0,1);
%
%X=zeros(n,leps,thismany);
%Y=X;
X=nan(leps,years,thismany);
Y=X;
thresh=.005;
yy=1:years;
for i=1:leps
    epsi=eps(i);
    for j=1:thismany
        [f,g,~]=finalSizeMulti(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,isdual,solvetype,numseed,epsi,1,tauend);
        %
        fx=f(:,burn+1:end); %g=g(burn+1:end);
        fx1=fx;
        %fsum=sum(fx,1);
        fsum=max(fx,[],1);
        fx(:,fsum<thresh)=[];
        yyi1=yy; yyi1(fsum<thresh)=[];
        lf=size(fx,2);
        %
        gx=g(:,burn+1:end); %g=g(burn+1:end);
        gx1=gx;
        %gsum=sum(gx,1);
        gsum=max(gx,[],1);
        gx(:,gsum<thresh)=[];
        yyi2=yy; yyi2(gsum<thresh)=[];
        lg=size(gx,2);
        %
        XX=zeros(n,lf);
        YY=zeros(n,lg);
        for k=1:lf
            %{
            xi=[(1:n)',fx(:,k)];
            xi=sortrows(xi,2);
            XX(:,k)=xi(:,1);
            %}
            cck=corrcoef(fx(:,k),fp);
            X(i,yyi1(k))=cck(2);
        end
        for k=1:lg
            %{
            yi=[(1:n)',gx(:,k)];
            yi=sortrows(yi,2);
            YY(:,k)=yi(:,1);
            %}
            cck=corrcoef(gx(:,k),fp);
            Y(i,yyi2(k))=cck(2);
        end
        %X(:,i,j)=var(XX,0,2);
        %Y(:,i,j)=var(YY,0,2);
    end
end
f=X;
g=Y;
save('MA1driftCC','X','Y')
end