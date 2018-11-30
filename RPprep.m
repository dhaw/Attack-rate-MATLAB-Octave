function [gamma,NN,n,K,L,Vi,Viover,kangle,h,w,rhohat,isflat,beta,boxLat,boxLong]=RPprep(lscan,R0,kangle,h,p)
%p is alpha
%Parameters:
aa=1;%Offset
gamma=1/10;
celldist=1;
w=kangle-h+1;

isflat=length(lscan);%Flat if just input density
if isflat==1
    sidea=5000;%>1 - if !d then must be column
    sideb=1;
    rho=h;
    lscan=rho*ones(sidea,sideb);
    rhohat=rho/w;
else
    [sidea,sideb]=size(lscan);
    rhohat=0;
end
%
%Population density:
[l1,l2]=size(lscan);
n=l1*l2;
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
%Vector of cell pops:
NN=lscan;
NN=reshape(NN,n,1);
NN=ceil(NN);
NN(NN<0)=NaN;
NN(isnan(NN)==1)=0;
%}
%%
%Fluscape:
K=1./(1+(r./aa).^(p));
%{
Nalpha=NN'.^alpha;
Njalpha=repmat(Nalpha,n,1);
Nione=repmat(NN,1,n);
K=K.*Njalpha.*Nione;
%}
%{
%Truscott:
K=1./(1+(r./13.5).^3.9)+.3*eye(n);
Nalpha=NN'.^.95;
Njalpha=repmat(Nalpha,n,1);
K=K.*Njalpha;
%}
K=K-diag(diag(K));%Zero diagonal elements
sumK=sum(K,2);
repK=repmat(sumK,1,n);
repK(repK==0)=1;
K=K./repK;
%
if isflat==1
    Vi=rho*ones(n,1);
else
    Vi=K*NN;
end
L=K+diag(Vi)/w;
%
%%
%R0 calculation:
%NGM with age only:
if isflat==1
    beta=R0*gamma*rhohat*kangle/(rho*(rhohat+1));
    Viover=0;
else
%Can't use quick method if pop not flat:
Sstart=repmat(NN,1,n);
Viover=repmat(Vi,1,n);
Viover(Viover==0)=1;
%
D=w/kangle*Sstart.*Viover.*L;
G=1/gamma*D;
d=eigs(G,1); R0a=max(d); beta=R0/R0a;
end
boxLat=sidea;
boxLong=sideb;
end