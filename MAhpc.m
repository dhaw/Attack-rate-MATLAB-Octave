function [f,g]=MAhpc%(C,Q,Qeven,fpand)
load('forMAhpc.mat')
%eps=(0:1:1);
%leps=length(eps);
chi=(0:1:1);
lchi=length(chi);
thismany=50;%Random ICs in fSMR
tauend=20;
burn=10;
%years=tauend-burn;
isdual=2;
solvetype=2;
numseed=10^(-8);
[gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,ages0]=prepFluAgeLocs(C,Qeven,0,1);
%X=zeros(n,thismany,leps);
%Y=zeros(thismany,leps);
Z1=zeros(n,thismany,lchi);
Z2=Z1;
%for i=1:leps
%epsi=eps(i);
epsi=.13;
for j=1:lchi
    chij=chi(j);
for k=1:thismany
    %FYI finalSizeMultiRank outputs changed
    [f,g]=finalSizeMulti2subRank(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,isdual,tauend,epsi,chij);
    %sample=(1:5:n);
    f=f(:,burn+1:end);
    fsum=sum(f,1);
    Z1(:,k,j)=var(f(:,fsum>.001),0,2);
    g=g(:,burn+1:end);
    gsum=sum(fg,1);
    Z2(:,k,j)=var(g(:,gsum>.001),0,2);
    %cc=corrcoef(f(:,end),fpand);
    %Y(j,i)=cc(2);
end
end
%end
f=Z1;
g=Z2;
save('MA2subZ','Z1','Z2')
%g=Y;
%}
%{
    %Bootstrap:
    A=datasample(s,f,200,2);%,'Replace',false)
    A=reshape(A,n,2,100);
    %Rank:
    X=[(1:n)',A];
    %}