function [f,g]=aggCells(D,k)
[a,b]=size(D);
ak=k*floor(a/k);
bk=k*floor(b/k);
a=ak; b=bk;
n2=zeros(a/k,b/k);
    for k1=1:k:a-k+1
        idxr=1+floor(k1/k);
        for k2=1:k:b-k+1
            idxc=1+floor(k2/k);
            NF=D(k1:k1+k-1,k2:k2+k-1);
            tempn=max(sum(NF(:)),1);% tempn(tempn==0)=1;
            n2(idxr,idxc)=tempn;
        end
    end
f=n2;