function [Df,xf,yf,Ds,Is]=fluscapeOnlyLocalAverage(lscan,fluscapeLocations)
[Df,xf,yf]=lscanplot(lscan,fluscapeLocations);

[a,b]=size(Df);
X=zeros(a+2,b+2);
X(2:end-1,2:end-1)=Df;
Y=Df+X(1:a,1:b)+X(1:a,2:b+1)+X(1:a,3:b+2);
Y=Y+X(2:a+1,1:b)+X(2:a+1,3:b+2);
Y=Y+X(3:a+2,1:b)+X(3:a+2,2:b+1)+X(3:a+2,3:b+2);
Y(2:end-1,2:end-1)=Y(2:end-1,2:end-1)/9;
Y(2:end-1,1)=Y(2:end-1,1)/6;
Y(2:end-1,end)=Y(2:end-1,end)/6;
Y(1,2:end-1)=Y(1,2:end-1)/6;
Y(end,2:end-1)=Y(end,2:end-1)/6;
Y(1,1)=Y(1,1)/4;
Y(1,end)=Y(1,end)/4;
Y(end,1)=Y(end,1)/4;
Y(end,end)=Y(end,end)/4;

[Ds,Is]=selectSubgrid(Y,xf,yf);
Ds=round(Ds);