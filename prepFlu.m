function [gamma,NN,n,nbar,na,NNbar,NNrep,maxN,cen,Kbar,K1,Cbar,betaS,betaI,betaD,ages0]=prepFlu(lscan,R0,stoch)%,U)
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
celldist=1;
%
Cnum=[6.92,.25,.77,.45;.19,3.51,.57,.2;.42,.38,1.4,.17;.36,.44,1.03,1.83];
Cdur=[3.88,.28,1.04,.49;.53,2.51,.75,.5;1.31,.8,1.14,.47;1,.85,.88,1.73];
C=Cnum.*Cdur;
%C=[6,.5;1,.5];%Steven's paper
Ca=C; Cb=C;
%Comment in to turn age off:
%Ca=ones(4); Cb=ones(4);
%
na=length(C);
a1=5;%1st age group - up to and including
a2=19;
a3=64;
%
%Population density:
[l1,l2]=size(lscan);
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
NN=lscan;
NN=reshape(NN,n,1);
NN=ceil(NN);
NN(NN<0)=NaN;
NN(isnan(NN)==1)=0;
%}
%NN=reshape(NNbar,n,4); NN=sum(NN,2);
%
%[maxN,cen]=max(NN);
%maxN=maxN(1);
%cen=cen(1);
%%
%Age:
%
v=[5,14,45,16]/80;
NNbar=[NN*v(1);NN*v(2);NN*v(3);NN*v(4)];
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
CC=[NN,cumsum(ones(n,1))];
CC=sortrows(CC,1);
CC=flipud(CC);
cen=CC(1,2);%1st arg-th most populous
maxN=CC(1,1);
%Redundant if input NNbar:
ages=sparse(n,maxN);
%
%No Urban/rural distinction:
%
K=1./(1+(r./aa).^(p));
Nalpha=NN'.^alpha;
Njalpha=repmat(Nalpha,n,1);
Nione=repmat(NN,1,n);
K=K.*Njalpha.*Nione;
sumK=sum(K,2);
repK=repmat(sumK,1,n);
repK(repK==0)=1;
K=K./repK;
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
Kbar(1:n,1:n)=eye(n); K1(1:n,1:n)=eye(n); 

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
%%
%beta_i:
%{
r=rand(n,1);
epsilon=.1; betai=(1-epsilon/2+epsilon*r); betai=repmat(betai,1,n);
%}
%ABM ages:
%{
amax=80;
for i=1:n
    Ni=NN(i);
    ages(i,1:Ni)=ceil((amax)*rand(1,Ni));
end
age1=sum(ages<=a1,2);
age2=sum(ages<=a2,2)-age1;
age3=sum(ages<=a3,2)-age1-age2;
age1=age1-sum(ages==0,2);
age4=NN-age1-age2-age3;
NNbar=[age1;age2;age3;age4];
%}