function [NN,r]=fluscapeNNr(Df,Ds)
celldist=1;
[l1,l2]=size(Ds);
n=l1*l2;
NN=reshape(Df,n,1);
ind=reshape(Ds,n,1);
NN(ind==0)=[];
ind=find(ind==1);

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
r=r(ind,ind);