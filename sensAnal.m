function [AS,AD,AI]=sensAnal(D)%[AS,AD,AI]
[l1,l2]=size(D); %n=l1*l2;
k=4;
ll=l1*l2-2*k*(l1+l2)+4*k^2;
R0=1.8;
param=(.1:.1:2); lp=length(param);
%
it=1;
AS=zeros(ll,lp); AI=AS; AD=AS;%,it
%B=zeros(n,lar,it);
%j=1;
for i=1:lp
    parami=param(i);
    [gamma,NN,n,nbar,na,NNbar,NNrep,maxN,cen,Kbar,K1,Cbar,betaS,betaI,betaD,ages0]=XprepFlu(D,R0,parami);
    %for j=1:it
    %
    try
        [zs,zs2]=ZfinalSizeAllMulti(gamma,n,nbar,na,NN,NNbar,NNrep,Kbar,K1,Cbar,betaS,0);
    catch ME
        zs=-1*ones(nbar,1);
    end
    zs=trimbyk(zs);
    AS(:,i)=zs;
    %}
    %
    try
        [zd,zd2]=ZfinalSizeAllMulti(gamma,n,nbar,na,NN,NNbar,NNrep,Kbar,K1,Cbar,betaD,1);
    catch ME
        zd=-1*ones(nbar,1);
    end
    zd=trimbyk(zd);
    AD(:,i)=zd;
    %}
    %
    count=1;
    while count>0 && count<10
    try
        [zi,zi2]=XfinalSizeAll(gamma,n,nbar,na,NN,NNbar,NNrep,Kbar,K1,Cbar,betaI,2);
        count=0;
    catch ME
        zi=-1*ones(n,1);
        count=count+1;
    end
    zi=trimbyk(zi);
    AI(:,i)=zi;
    %}
    end
end