function f=finalSize1YDE(gamma,n,nbar,na,NN,NNbar,NNrep,Kout,Kin,K1,Cbar,beta)%1 year, delta/epsilon
NN0=NNrep; NN0(NNrep==0)=1;
%
Mj=Kout'*NNbar;%Approximation
Mjover=1./Mj; Mjover(Mj==0)=0; Mjover=repmat(Mjover',nbar,1);
%
L=beta*(K1.*Mjover)*(Kin'.*Cbar);
zn=zeros(nbar,1);
Z0=zn;
f=XODEsolveAllMulti(gamma,NN,n,nbar,na,NNbar,NNrep,NN0,Kout,Kin,K1,Cbar,zn,Z0,beta,Mjover,L);
end

function f=XODEsolveAllMulti(gamma,NN,n,nbar,na,NNbar,NNrep,NN0,Kout,Kin,K1,Cbar,zn,Z0,beta,Mjover,L)
%nbar=length(NNbar);
t0=0; tend=360;
%
NNprob=NNbar/sum(NN);
numseed=10^(-8); seed=numseed*NNprob;
icR=Z0.*NNrep;
y0=[NNbar-icR;zn;icR];

[tout,yout]=ode45(@(t,y)integr8all(t,y,beta,gamma,n,nbar,NNbar,NN0,Kout,Kin,K1,Cbar,seed,Mjover,L),[t0,tend],y0);
rec0=yout(end,2*nbar+1:end);
Zsol=rec0'./NN0;
f=sum(reshape(Zsol,n,na),2);
end

function f=integr8all(t,y,beta,gamma,n,nbar,NNbar,NN0,Kout,Kin,K1,Cbar,seed,Mjover,L)
S=y(1:nbar);
I=y(nbar+1:2*nbar);
foi=L*I+seed;
Sdot=-S.*foi;
Idot=S.*foi-gamma*I;
Rdot=gamma*I;
f=[Sdot;Idot;Rdot];
end