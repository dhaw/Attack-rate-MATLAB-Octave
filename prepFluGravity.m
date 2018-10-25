function [gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,ages0]=prepFluGravity(lscan,R0,stoch,K,NN)%,U)
%stoch=1 for SCM(/ABM)
%Parameters:
gamma=1/2.6;
a1immobile=0;
normaliseKernel=1;
%%
%
Cnum=[6.92,.25,.77,.45;.19,3.51,.57,.2;.42,.38,1.4,.17;.36,.44,1.03,1.83];
Cdur=[3.88,.28,1.04,.49;.53,2.51,.75,.5;1.31,.8,1.14,.47;1,.85,.88,1.73];
C=Cnum;%.*Cdur;
C=ones(4);
%}
%Comment in to turn age off:
%C=ones(4);
%Ca=C; Cb=C;
%If plotting curves, can trimbyk before finding max/min **1**
na=length(C);
a1=5;%1st age group - up to and including
a2=19;
a3=64;
%%
%Population density:
n=length(NN);
nbar=n*na;
%%
%Age:
v=[5,14,45,16]/80;
NNbar=[NN*v(1);NN*v(2);NN*v(3);NN*v(4)];
if stoch==1
    NNbar=round(NNbar);
    NN=sum(reshape(NNbar,n,na),2);
end
%NNbar=[.2*NN;.8*NN];
NNrep=repmat(NN,na,1);
Sstart=repmat(NNbar,1,nbar);
SstartFrac=NNbar./NNrep; SstartFrac(NNrep==0)=0; SstartFrac=repmat(SstartFrac,1,nbar);
%ABM: ages in here (code at bottom)
%%
%CC=[trimbyk(NN),trimbyk(cumsum(ones(n,1)))];   **1**
CC=[NN,cumsum(ones(n,1))];
CC(CC(:,1)==0,:)=[];
CC=sortrows(CC,1);
CC=flipud(CC);
%cen=CC(1,2);%1st arg-th most populous
minNind=CC(end,2);
maxNind=CC(1,2);
maxN=CC(1,1);
%Redundant if input NNbar:
ages=sparse(n,maxN);
sumK=sum(K,2);
repK=repmat(sumK,1,n);
repK(repK==0)=1;
if normaliseKernel==1
    K=K./repK;
end
%
%
%Expand K and C:
Kbar=repmat(K,na,na);
keye=eye(n);
kxeye=ones(n)-keye;
Cbar=kron(Ca,keye)+kron(Cb,kxeye);
%%
%R0 calculation:
K1=kron(eye(na),K);
%Children immobile:
if a1immobile==1
    Kbar(1:n,1:n)=eye(n); K1(1:n,1:n)=eye(n);
end
%NGM with age only:
if normaliseKernel==1
    Ntot=sum(NN);
    Asum=reshape(NNbar,n,na); Asum=sum(Asum,1);
    X=repmat(Asum',1,na).*C/Ntot/gamma;
    d=eigs(X,1); R0a=max(d); betaS=R0/R0a; betaD=betaS; betaI=betaS;
else
%Can't use quick method if K not normalised:
Mj=Kbar'*NNbar;
Mj(Mj==0)=1;
Mjover=1./Mj;
Mjover=repmat(Mjover',nbar,1);

DS=Sstart.*Kbar.*Mjover.*Cbar;
GS=1/gamma*DS;

DI=SstartFrac.*Kbar'.*Cbar;
GI=1/gamma*DI;

DD=(Sstart.*K1.*Mjover)*(Kbar'.*Cbar);
GD=1/gamma*DD;

d=eigs(GS,1); R0a=max(d); betaS=R0/R0a;
d=eigs(GI,1); R0a=max(d); betaI=R0/R0a;
d=eigs(GD,1); R0a=max(d); betaD=R0/R0a;
end
ages0=ages;
end