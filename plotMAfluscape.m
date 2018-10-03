function f=plotMAfluscape(Asus,Arec,NN,NNbar,Ds,Is)
%Outputs "best guess" of mean H3 attack rate
n=length(Is);
na=4;
nbar=na*n;
%locs=reshape(Ds,n,1);
%inds=reshape(Is,n,1);
%inds(locs==0)=[];
inds=Is;
locind=find(inds>0);
locno=inds(locind);
numlocs=length(locno);

v1=[locind;locind+n;locind+2*n;locind+3*n];
v2=v1+nbar;
v3=v2+nbar;
%v4=v3+nbar;
Asus=Asus([v1;v2;v3],:);
Arec=Arec([v1;v2],:);
NN=NN(locind,:);
NNbar=NNbar(v1,:);
plotMAspace(Asus,Arec,NN,NNbar)
X=Arec(na*numlocs+1:end,:);%H3 only

Xsum=sum(X,1);

X=nanmean(X(:,Xsum>0),2);
X=reshape(X,numlocs,na); X=sum(X,2)./NN; X(NN==0)=0;
f=[locno,X];