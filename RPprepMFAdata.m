function [gamma,NN,n,K,L,Vi,Viover,kangle,h,w,rhohat,isflat,beta,boxLat,boxLong]=RPprepMFAdata(X,R0,kangle,h,p)
%HH,rep,lat,long
factor=;%Long/lat to km
lh=size(X,1);
[x1,x2]=meshgrid(X(:,1),X(:,1));
x=(x1-x2).^2;
[y1,y2]=meshgrid(X(:,2),X(:,2));
y=(y1-y2).^2;
r=sqrt(x+y)*factor;
%Coarse grain for boxLat/boxLong? - for plot