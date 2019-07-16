function [f,g,h]=MAhpc1rank%(C,Q,Qeven,fpand)
load('forMAhpc.mat')
eps=(0:.01:1);
leps=length(eps);
%chi=(0:.01:1);
thresh=10^(-4);
%
thismany=1;%Random ICs in fSMR - loop length
tauend=100;
burn=50;
%years=tauend-burn;
isdual=3;
solvetype=2;
numseed=10^(-8);
[gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,ages0]=prepFluAgeLocs(C,Qeven,0,1);
X=zeros(n,tauend-burn,thismany);
Y=zeros(n,thismany);
Z=zeros(tauend-burn,thismany);
%for i=1:leps
%epsi=eps(i);
for j=1:thismany
    %try
    [f,g]=finalSizeMulti(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,isdual,solvetype,numseed);
    gx=g(:,burn+1:end); %g=g(burn+1:end);
    gx1=gx;
    gsum=sum(gx,1);
    inds=find(gsum<thresh);
    %gx(:,inds)=[];
    %lg=size(gx,2);
    
    K=K1(1:n,1:n);
    [v,lam]=eig(K);
    lam=diag(lam); [~,lam1ind]=max(lam);
    v1=v(:,lam1ind);
    
    for i=1:tauend-burn
        if gsum(i)<thresh
            Z(i,j)=0;
        else
        xi=[(1:n)',gx(:,i)];
        xi=sortrows(xi,2);
        X(:,i,j)=xi(:,2);
        
        zi=corrcoef(gx1(:,i),v1);
        Z(i,j)=zi(1,2);
        end
    end
    Y(:,j)=var(X(:,:,j),0,2);
    %catch ME
    %end
end
%end
Xranks=X;
Yvars=Y;
Zccs=Z;
f=Xranks;
g=Yvars;
h=Zccs;
%save('MAtest','Xranks','Yvars','Zccs')