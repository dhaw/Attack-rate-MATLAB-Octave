function f=fluscapeSeed(Ds)
[a,b]=size(Ds);
n=a*b;
S=zeros(a,b);
S(11:15,end-14:end-10)=ones(5);
S=S.*Ds;
sample=reshape(S,n,1);
Nind=reshape(Ds,n,1);
sample(Nind==0)=[];
f=sample;