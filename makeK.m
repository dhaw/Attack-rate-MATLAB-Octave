function [K1,K2,K3,K4,NN]=makeK(lscan,r)
lscan=round(lscan);
lscan(lscan<0)=0;
%Kernel=1(G) 2(OG) 3(R) 4(OR)
celldist=1;%km
alpha1=0.514160073342821;
p1=2.66171686998094;
%
aa2=10^(-0.231446830279379);
alpha2=0.518017124824123;
p2=2.74310227108717;
%
aa4=3.17785167961894;
%%
NN=lscan;
NN(isnan(NN)==1)=0;
[l1,l2]=size(NN);
n=l1*l2;
%Comment out if input population vector and distances:
if r==1
%
NN=reshape(NN,n,1);
NN=ceil(NN);
NN(NN<0)=NaN;
%
%Kernel - distances only:
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
end
%}
%%
K1=1./(1+r.^p1);
Nalpha=NN'.^alpha1;
Njalpha=repmat(Nalpha,n,1);
%Nione=repmat(NN,1,n); Nione ignorable if going to row normalise
K1=K1.*Njalpha;%.*Nione;
%%
K2=1./(1+(r./aa2).^p2);
Nalpha=NN'.^alpha2;
Njalpha=repmat(Nalpha,n,1);
%Nione=repmat(NN,1,n);
K2=K2.*Njalpha;%.*Nione;
%%
S=zeros(n); Sp=S;
%nvec=(1:n);
Srow=zeros(1,n); Sprow=Srow;
    for i=1:n
        x=zeros(n,1);
        %
        rowi=r(i,:);
        [b,ind]=sort(rowi);%b=rowi(ind);
        %backind=nvec(ind);%Indices to shuffle back to original order - rowi=b(backind)
        Norder=NN; Norder(i)=0;
        Norder=Norder(ind);%Cell pops in order of distance
        [C,ib,ic]=unique(b);%C=b(ib), b=C(ic)
        %For each new C, which cells have distance less than
        upto=0;
        for j=2:length(C)
            findj=find(ic==j);%ic is ordered, so there may be a quicker way
            upto=upto+sum(Norder(findj));
            x(findj)=upto;
        end
        %For offset:
        off=find(C<=aa4);%All distances within offset
        off=off(end);%Index in C of furthest distance within offset
        indoff=find(ic<=off);%indices in b/Nb/x (distances) of cells within offset
        indoff(indoff==1)=[];
        xoff=x(indoff(end));%x-vale of furthest cell within offset %ic(off) C(k)=b(ib(k))
        y=x; y(indoff)=xoff;%Replace x-values within offset with xoff
        x=x-Norder; y=y-Norder;%Subtract Nj from each term: %x(i)=0;
        Srow(ind)=x;%Re-shuffle x to original cell order
        Sprow(ind)=y;
        S(i,:)=Srow'; %S(i,i)=0; 
        Sp(i,:)=Sprow'; %Sp(i,i)=0;
        %}
    end
Ni=repmat(NN,1,n); Nj=Ni'; NiNj=Ni.*Nj;
denom3=(S+Ni).*(S+Ni+Nj);
K3=NiNj./denom3; K3(denom3==0)=0;
%%
denom4=(Sp+Ni).*(Sp+Ni+Nj);
K4=NiNj./denom4; K4(denom4==0)=0;
%{
        rowi=r(i,:);
        [C,ia,ic]=unique(rowi); %C=rowi(ia), rowi=C(ic)
        %For each new C, which cells have distance less than
        [Corder,ind]=sort(C);%Corder=C(ind);
        ia=ia(ind); ic=ic(backind);
        
        Nb=Corder(ic);
        upto=0;
        for j=2:lc
            findj=find(newic==j);%Find cells up to j'th distance away
            upto=upto+sum(Nb(findj));%Add populations
            x(findj)=upto;%put into ranked S_ij vector
        end
        %For offset:
        off=find(C>=aa4); off=off(1); xoff=x(ic(off));
        y=x; y(y<xoff)=xoff;
        %
        x=x-Nb; y=y-Nb;%Subtract Nj from each entry %x(i)=0;
        S(i,:)=x(ind); S(i,i)=0;%Unrank
        Sp(i,:)=y(ind); Sp(i,i)=0;
        %}