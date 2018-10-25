function [gamma,n,betaS,betaI,betaD,K]=prepFluGravityNoAge(R0,K,NN)%,U)
%NN must be rounded
%Parameters:
gamma=1/2.6;
normaliseKernel=1;
%%
n=length(NN);
Sstart=repmat(NN,1,n);
SstartFrac=NN; SstartFrac(NN==0)=0; SstartFrac=repmat(SstartFrac,1,n);
%ABM: ages in here (code at bottom)
%%
sumK=sum(K,2); sumK(sumK==0)=1;
repK=repmat(sumK,1,n);
repK(repK==0)=1;
if normaliseKernel==1
    K=K./repK;
end
%NGM with age only:
if normaliseKernel==1
    betaS=R0*gamma;
    betaI=betaS;
    betaD=betaS;
else
%Can't use quick method if K not normalised:
Mj=Kbar'*NN;
Mj(Mj==0)=1;
Mjover=1./Mj;
Mjover=repmat(Mjover',n,1);

DS=Sstart.*K.*Mjover;
GS=1/gamma*DS;

DI=SstartFrac.*K';
GI=1/gamma*DI;

DD=(Sstart.*K.*Mjover)*K';
GD=1/gamma*DD;

d=eigs(GS,1); R0a=max(d); betaS=R0/R0a;
d=eigs(GI,1); R0a=max(d); betaI=R0/R0a;
d=eigs(GD,1); R0a=max(d); betaD=R0/R0a;
end
end