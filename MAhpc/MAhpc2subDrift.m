function [f,g,h]=MAhpc2subDrift%(lscan,fluscapeLocations)
thismany=50;%Random ICs in fSMR - loop length
tauend=1000;
burn=500;
years=tauend-burn;
isdual=0;
%solvetype=2;
numseed=10^(-5);
eps=(0:.1:1);
leps=length(eps);
chi=.2;
filename='MA1driftCCchip2';
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
[fp,gp,~]=finalSizeMulti(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,isdual,2,numseed,0,0,1);
repNN=repmat(NN,1,years);
%
%X=zeros(n,leps,thismany);
%Y=X;
X=nan(leps,years,thismany);
Y=X;
Z=X;
thresh=.005;
yy=1:years;
for i=1:leps
    epsi=eps(i);
    for j=1:thismany
        [f,g]=finalSizeMulti2sub(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,isdual,numseed,epsi,chi,1,tauend);
        %
        R1=g(1:nbar,burn+1:end);
        R1=R1(1:n,:)+R1(n+1:2*n,:)+R1(2*n+1:3*n,:)+R1(3*n+1:end,:);
        fx=R1./repNN; R1(repNN==0)=0;
        R2=g(nbar+1:end,burn+1:end);
        R2=R2(1:n,:)+R2(n+1:2*n,:)+R2(2*n+1:3*n,:)+R2(3*n+1:end,:);
        gx=R2./repNN; R2(repNN==0)=0;
        hx=fx+gx;
        
        fsum=max(fx,[],1);
        fx1=fx;
        fx(:,fsum<thresh)=[];
        yyi1=yy; yyi1(fsum<thresh)=[];
        lf=size(fx,2);
        %
        gsum=max(gx,[],1);
        gx1=gx;
        gx(:,gsum<thresh)=[];
        yyi2=yy; yyi2(gsum<thresh)=[];
        lg=size(gx,2);
        %
        hsum=max(hx,[],1);
        hx1=hx;
        hx(:,hsum<thresh)=[];
        yyi3=yy; yyi3(hsum<thresh)=[];
        lh=size(hx,2);
        for k=1:lf
            %{
            xi=[(1:n)',fx(:,k)];
            xi=sortrows(xi,2);
            XX(:,k)=xi(:,1);
            %}
            cck=corrcoef(fx(:,k),fp);
            X(i,yyi1(k),j)=cck(2);
        end
        for k=1:lg
            %{
            yi=[(1:n)',gx(:,k)];
            yi=sortrows(yi,2);
            YY(:,k)=yi(:,1);
            %}
            cck=corrcoef(gx(:,k),fp);
            Y(i,yyi2(k),j)=cck(2);
        end
        for k=1:lh
            %{
            yi=[(1:n)',gx(:,k)];
            yi=sortrows(yi,2);
            YY(:,k)=yi(:,1);
            %}
            cck=corrcoef(hx(:,k),fp);
            Z(i,yyi3(k),j)=cck(2);
        end
        %X(:,i,j)=var(XX,0,2);
        %Y(:,i,j)=var(YY,0,2);
        %h1prop=fx1./(fx1+gx1);
        %Z(i,:,j)=nanvar(h1prop,[],2);
    end
end
f=X;
g=Y;
h=Z;
save(filename,'X','Y','Z')
end