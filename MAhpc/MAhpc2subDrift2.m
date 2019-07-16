function [f,g,h]=MAhpc2subDrift2%(lscan,fluscapeLocations)
thismany=2;%Random ICs in fSMR - loop length
tauend=20;
burn=10;
years=tauend-burn;
isdual=0;
%solvetype=2;
numseed=10^(-5);
eps=(0:.1:.2);
del=(-.1:.1:.1);
%
leps=length(eps);
ldel=length(del);
chi=.2;%Change filename also********
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
xnan=nan(thismany,years);%(leps,years,thismany);
X=cell(leps,ldel);
[X{:,:}]=deal(xnan);
xnanRow=nan(1,years);
%X=nan(leps,ldel,thismany,years);
Y=X;
Z=X;
thresh=.005;
yy=1:years;
parfor i1=1:leps
epsi1=eps(i1);
for i2=1:ldel
deli2=del(i2);
    for j=1:thismany
        x=xnanRow; y=xnanRow; z=xnanRow;
        [f,g]=finalSizeMulti2sub(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,isdual,numseed,epsi1,deli2,chi,1,tauend);
        %
        R1=g(1:nbar,burn+1:end);
        R1=R1(1:n,:)+R1(n+1:2*n,:)+R1(2*n+1:3*n,:)+R1(3*n+1:end,:);
        fx=R1./repNN; %R1(repNN==0)=0;
        R2=g(nbar+1:end,burn+1:end);
        R2=R2(1:n,:)+R2(n+1:2*n,:)+R2(2*n+1:3*n,:)+R2(3*n+1:end,:);
        gx=R2./repNN; %R2(repNN==0)=0;
        hx=fx+gx;
        
        fsum=max(fx,[],1);
        %fx1=fx;
        fx(:,fsum<thresh)=[];
        yyi1=yy; yyi1(fsum<thresh)=[];
        lf=size(fx,2);
        %
        gsum=max(gx,[],1);
        %gx1=gx;
        gx(:,gsum<thresh)=[];
        yyi2=yy; yyi2(gsum<thresh)=[];
        lg=size(gx,2);
        %
        hsum=max(hx,[],1);
        %hx1=hx;
        hx(:,hsum<thresh)=[];
        yyi3=yy; yyi3(hsum<thresh)=[];
        lh=size(hx,2);
        for k=1:lf
            cck=corrcoef(fx(:,k),fp);
            x(1,yyi1(k))=cck(2);%x(j,
            X{i1,i2}(j,:)=x;
            %X(i1,i2,j,yyi1(k))=cck(2);
        end
        for k=1:lg
            cck=corrcoef(gx(:,k),fp);
            y(1,yyi2(k))=cck(2);
            Y{i1,i2}(j,:)=y;
            %Y(i1,i2,j,yyi2(k))=cck(2);
        end
        for k=1:lh
            cck=corrcoef(hx(:,k),fp);
            z(1,yyi3(k))=cck(2);
            Z{i1,i2}(j,:)=z;
            %Z(i1,i2,j,yyi3(k))=cck(2);
        end
    end
%X{i1,i2};
%Y{i1,i2}=y;
%Z{i1,i2}=z;
%save(filename,'X','Y','Z')
end
end
f=X;
g=Y;
h=Z;
save(filename,'X','Y','Z')
end