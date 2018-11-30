function f=MAhpc2%(C,Q,Qeven,fpand)
load('forMAhpc.mat')
eps=(0:.01:1);
leps=length(eps);
%chi=(0:.01:1);
%
thismany=5;%Random ICs in fSMR
tauend=100;
burn=50;
%years=tauend-burn;
isdual=3;
solvetype=2;
numseed=10^(-8);
[gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,ages0]=prepFluAgeLocs(C,Qeven,0,1);
X=zeros(n,tauend-burn,thismany);
Y=zeros(n,thismany);
%for i=1:leps
%epsi=eps(i);
for j=1:thismany
    %try
    [f,g]=finalSizeMulti(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,isdual,solvetype,numseed);
    gx=g(:,burn+1:end); %g=g(burn+1:end);
    for i=1:tauend-burn
        xi=[(1:n)',gx(:,i)];
        xi=sortrows(xi,2);
        X(:,i)=xi(:,2);
    end
    Y(:,j)=var(X(:,:,j),0,2);
    %catch ME
    %end
end
%end
f=X;
g=Y;
save('MAtest','X','Y')