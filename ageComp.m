<<<<<<< HEAD
%function f=ageComp(B1,B2,B3,B4)
function f=ageComp(C1,C2,C3,C4,gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,ages0)
%{
C1=aggCells(B1,2);
C2=aggCells(B2,2);
C3=aggCells(B3,2);
C4=aggCells(B4,2);
[a,b]=size(C1);
[gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,ages0]=prepFluAge(C1,C2,C3,C4,1.8,0);
%}
[fd,gd]=finalSizeMultiAge(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaD,1,2,10^(-5));
%[fs,gs]=finalSizeMultiAge(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,0,2,10^(-5));
%[fi,gi]=finalSizeMultiAge(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaI,2,2,10^(-5));

n1=reshape(C1,a*b,1); %n1=trimbyk(n1,a,b);
n2=reshape(C2,a*b,1); %n2=trimbyk(n2,a,b);
n3=reshape(C3,a*b,1); %n3=trimbyk(n3,a,b);
n4=reshape(C4,a*b,1); %n4=trimbyk(n4,a,b);
%{
Fd=trimbyk(fd,a,b); Gd=trimbyk(gd,a,b);
Fs=trimbyk(fs,a,b); Gs=trimbyk(gs,a,b);
Fi=trimbyk(fi,a,b); Gi=trimbyk(gi,a,b);
%}
plotAR(g(:,1),n1,a,b,'0-4')
plotAR(g(:,1),n1,a,b,'0-4')
plotAR(g(:,1),n1,a,b,'0-4')
plotAR(g(:,1),n1,a,b,'0-4')
%plotAR(fd,NN,a,b,'D-mobility')
%plotAR(fs,NN,a,b,'S-mobility')
=======
%function f=ageComp(B1,B2,B3,B4)
function f=ageComp(C1,C2,C3,C4,gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,ages0)
%{
C1=aggCells(B1,2);
C2=aggCells(B2,2);
C3=aggCells(B3,2);
C4=aggCells(B4,2);
[a,b]=size(C1);
[gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,ages0]=prepFluAge(C1,C2,C3,C4,1.8,0);
%}
[fd,gd]=finalSizeMultiAge(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaD,1,2,10^(-5));
%[fs,gs]=finalSizeMultiAge(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,0,2,10^(-5));
%[fi,gi]=finalSizeMultiAge(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaI,2,2,10^(-5));

n1=reshape(C1,a*b,1); %n1=trimbyk(n1,a,b);
n2=reshape(C2,a*b,1); %n2=trimbyk(n2,a,b);
n3=reshape(C3,a*b,1); %n3=trimbyk(n3,a,b);
n4=reshape(C4,a*b,1); %n4=trimbyk(n4,a,b);
%{
Fd=trimbyk(fd,a,b); Gd=trimbyk(gd,a,b);
Fs=trimbyk(fs,a,b); Gs=trimbyk(gs,a,b);
Fi=trimbyk(fi,a,b); Gi=trimbyk(gi,a,b);
%}
plotAR(g(:,1),n1,a,b,'0-4')
plotAR(g(:,1),n1,a,b,'0-4')
plotAR(g(:,1),n1,a,b,'0-4')
plotAR(g(:,1),n1,a,b,'0-4')
%plotAR(fd,NN,a,b,'D-mobility')
%plotAR(fs,NN,a,b,'S-mobility')
>>>>>>> acfad3e4da850e525286e08138beef9f06a8c2cc
%plotAR(fi,NN,a,b,'I-mobility')