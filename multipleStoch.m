function [f,g]=multipleStoch(D,R0)%(A,B,NN);
[a,b]=size(D);
ncut=(a-8)*(b-8);
runs=1;
A=zeros(ncut,runs);
B=A;
[gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,ages0]=prepFlu(D,R0,1);
for i=1:runs
    [fs,x]=finalSizeMulti(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,0,3,10^(-4));
    [fi,x]=finalSizeMulti(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaI,2,3,10^(-4));
    A(:,i)=trimbyk(fs,a,b);
    B(:,i)=trimbyk(fi,a,b);
end
f=A;
g=B;
%{
NN=trimbyk(NN);
nfind=find(NN==0);
A(nfind,:)=[]; B(nfind,:)=[]; NN(nfind)=[];
Amean=mean(A,2); Apct=prctile(A,[.23,.27],2);
Bvar=mean(B,2); ABpct=prctile(B,[.23,.27],2);
figure
%}