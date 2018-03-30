<<<<<<< HEAD
function f=trimbyk(A,a,b)
k=4;
%[a,b]=size(A);
%a=33;
%b=55;

ind=(1:a*b)';
ind=reshape(ind,a,b);
ind=ind(k+1:end-k,k+1:end-k);
ind=reshape(ind,(a-2*k)*(b-2*k),1);
f=A(ind,:,:);
%A=reshape(A,a,b);
%A=A(k+1:end-k,k+1:end-k);
%A=reshape(A,(a-2*k)*(b-2*k),1);
=======
function f=trimbyk(A,a,b)
k=4;
%[a,b]=size(A);
%a=33;
%b=55;

ind=(1:a*b)';
ind=reshape(ind,a,b);
ind=ind(k+1:end-k,k+1:end-k);
ind=reshape(ind,(a-2*k)*(b-2*k),1);
f=A(ind,:,:);
%A=reshape(A,a,b);
%A=A(k+1:end-k,k+1:end-k);
%A=reshape(A,(a-2*k)*(b-2*k),1);
>>>>>>> acfad3e4da850e525286e08138beef9f06a8c2cc
%f=A;