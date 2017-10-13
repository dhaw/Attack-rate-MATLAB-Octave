function [gamma,NN,n,nbar,na,NNbar,NNrep,Kout,Kin,K1,Cbar,beta]=prepFluDE(lscan,R0,delta,eps)%Delta Epsilon
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
%C=[6,.5;1,.5];
%Comment in to turn age off:
%C=ones(4);
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
%%
CC=[NN,cumsum(ones(n,1))];
CC=sortrows(CC,1);
CC=flipud(CC);
cen=CC(1,2);%1st arg-th most populous
maxN=CC(1,1);
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
%Delta epsilon:
E=eye(n);
A=(1-delta)*E+delta*K;
B=(1-eps)*E+eps*K;
%
%Expand K and C:
Abar=repmat(A,na,na); K1=kron(eye(na),A);
Bbar=repmat(B,na,na);
Cbar=kron(C,ones(n));
Kout=Abar; Kin=Bbar;
%%
%Age:
%Redundant if input NNbar:
ages=sparse(n,maxN);
%
v=[5,14,45,16]/80;
NNbar=[NN*v(1);NN*v(2);NN*v(3);NN*v(4)];
%NNbar=round(NNbar);
%NNbar=[.2*NN;.8*NN];
NNrep=repmat(NN,na,1);
%beta=R0*gamma;
%
Sstart=repmat(NNbar,1,nbar);
SstartFrac=NNbar./NNrep; SstartFrac(NNrep==0)=0; SstartFrac=repmat(SstartFrac,1,nbar);
%%
Mj=Kout'*NNbar; Mj=sum(Mj,2); Mjover=1./Mj; Mjover(Mj==0)=1; Mjover=repmat(Mjover',nbar,1);
foi=(K1.*Mjover)*(Kin'.*Cbar);
GD=1/gamma*Sstart.*foi;
d=eigs(GD,1); R0a=max(d); beta=R0/R0a;
%}