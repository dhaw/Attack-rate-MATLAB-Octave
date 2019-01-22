function f=prepRasterCsv(B)
B=B(:,2:4);
B(isnan(B(:,1))==1,:)=[];
%If 100m x 100m:
%
B(:,1)=(B(:,1)-min(B(:,1)))/100+1;
B(:,2)=(B(:,2)-min(B(:,2)))/100+1;
%}
%{
B(:,1)=(B(:,1)-min(B(:,1)))/1000+1;
B(:,2)=(B(:,2)-min(B(:,2)))/1000+1;
%}
Bsp=sparse(sparse(B(:,2),B(:,1),B(:,3)));%round
Bfull=full(Bsp);
f=flipud(Bfull);