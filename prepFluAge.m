function [gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,ages0]=prepFluAge(D1,D2,D3,D4,R0,stoch,delta)%,U)
%stoch=1 for SCM(/ABM)
%Parameters:
aa=.58;
aaR=.91;
aaU=0;
alpha=.52;
alphaR=.56;
alphaU=.4;
p=2.72;
pR=2.84;
pU=2.7;
gamma=1/2.6;
celldist=.2;%km
a1immobile=1;
normaliseKernel=1
%{
Cnum=[6.92,.25,.77,.45;.19,3.51,.57,.2;.42,.38,1.4,.17;.36,.44,1.03,1.83];
Cdur=[3.88,.28,1.04,.49;.53,2.51,.75,.5;1.31,.8,1.14,.47;1,.85,.88,1.73];
C=Cnum.*Cdur;
%}
%UK (from vaxedemic):
C=[37.4622640266199,13.2337799407673,9.35866526693108,5.27807222067513;17.2304141889828,98.1983003738366,17.0186152145963,10.1131975048866;9.46784315245156,9.4416088929148,16.22285757548,5.7675253611147;1.38284918679668,1.26680284573205,1.08367881504336,3.88324564380799];

%C=[6,.5;1,.5];%Steven's paper
Ca=C; Cb=C;
%Comment in to turn age off:
%Ca=ones(4); Cb=ones(4);
%If plotting curves, can trimbyk before finding max/min **1**
%
na=length(C);
%
%Population density:
[l1,l2]=size(D1);
n=l1*l2;
nbar=n*na;
%beta_i in here (code at bottom)
%%
%Kernel:
L=cumsum(ones(n,1));
L=reshape(L,[l1,l2]);%l1 column, labelled downwards
L=sparse(L);
[L1,L2,L3]=find(L);
L=[L1,L2,L3];
L=sortrows(L,3);
[x1,x2]=meshgrid(L(:,1),L(:,1));
x=(x1-x2).^2;
[y1,y2]=meshgrid(L(:,2),L(:,2));
y=(y1-y2).^2;
r=sqrt(x+y)*celldist;
%%
%
NN1=reshape(D1,n,1); NN2=reshape(D2,n,1); NN3=reshape(D3,n,1); NN4=reshape(D4,n,1);
NN=sum([NN1,NN2,NN3,NN4],2);
NN=ceil(NN);
NN(NN<0)=NaN;
NN(isnan(NN)==1)=0;
%}
%%
%Age:
%
NNbar=[NN1;NN2;NN3;NN4];
if stoch==1
    NNbar=round(NNbar);
end
NN=reshape(NNbar,n,4); NN=sum(NN,2);%Added********
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
minNind=CC(end,2);
maxNind=CC(1,2);
maxN=CC(1,1);
%Redundant if input NNbar:
ages=0;%sparse(n,maxN);
%
%No Urban/rural distinction:
%
K=1./(1+(r./aa).^(p));
Nalpha=NN'.^alpha;
Njalpha=repmat(Nalpha,n,1);
Nione=repmat(NN,1,n);
K=K.*Njalpha.*Nione;
%}
%
%Truscott:
K=1./(1+(r./13.5).^3.9)+.3*eye(n);
Nalpha=NN'.^.95;
Njalpha=repmat(Nalpha,n,1);
K=K.*Njalpha;
%}
%
sumK=sum(K,2);
repK=repmat(sumK,1,n);
repK(repK==0)=1;
if normaliseKernel==1
    K=K./repK;
end
%}
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
    %Kbar(1:n,1:n)=eye(n); K1(1:n,1:n)=eye(n);
    Kdelta=delta*K+(1-delta)*eye(n);
    Kbar(n+1:2*n,n+1:2*n)=Kdelta; K1(n+1:2*n,n+1:2*n)=Kdelta;
    Kbar(3*n+1:4*n,3*n+1:4*n)=Kdelta; K1(3*n+1:4*n,3*n+1:4*n)=Kdelta;
end
%NGM with age only - for I-mobility with approximation:
%{
Ntot=sum(NN);
Asum=reshape(NNbar,n,na); Asum=sum(Asum,1);
X=repmat(Asum',1,na).*C/Ntot/gamma;
d=eigs(X,1); R0a=max(d); betaS=R0/R0a; betaD=betaS; betaI=betaS;
%}
%
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
%}
ages0=ages;
end