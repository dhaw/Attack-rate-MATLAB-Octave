function [f,g,h]=MAhpcDelay
load('forMAhpc.mat')
delay=(-5:5);
ld=length(delay);
%
thismany=10;%Random ICs in fSMR
numseed=8;
eps=0;
del=0;
cross=1;
randic=0;
tauend=100;
burn=20;
%years=tauend-burn;
isdual=2;
solvetype=2;
[gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,ages0]=prepFluAgeLocs(C,Qeven,0,1);
%[fp,gp,~]=finalSizeMulti(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,isdual,2,numseed,0,0,1);
NN0=NN; NN0(NN==0)=1; NN0=repmat(NN0,1,tauend-burn);
%
F=nan(n,tauend-burn,thismany);
G=F;
bollocks=zeros(thismany,1);
for i=1:thismany;%ld
    %delayi=delay(i);
    delayi=5;%Max deviation from numseed
    %for j=1:thismany
    %try
        [f,g]=finalSizeMulti2subDelay(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,isdual,numseed,eps,del,cross,randic,tauend,delayi);
        %plotMA2sub(fi,gi,NN,NNbar,1,20);
        %
        R1=g(1:nbar,burn+1:end);
        R1=R1(1:n,:)+R1(n+1:2*n,:)+R1(2*n+1:3*n,:)+R1(3*n+1:end,:);
        fx=R1./NN0; %R1(repNN==0)=0;
        R2=g(nbar+1:end,burn+1:end);
        R2=R2(1:n,:)+R2(n+1:2*n,:)+R2(2*n+1:3*n,:)+R2(3*n+1:end,:);
        gx=R2./NN0; %R2(repNN==0)=0;
        %hx=fx+gx;
        F(:,:,i)=fx;
        G(:,:,i)=gx;
    %catch ME
        %bollocks(i)=1;
    %end
        %}
    %end
end
f=F;
g=G;
h=bollocks;
save('MAhpcDelaymax8.mat','F','G','bollocks')