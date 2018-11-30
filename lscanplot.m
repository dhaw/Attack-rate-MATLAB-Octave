function [Z2,x,y]=lscanplot(lscan,locs)
truncateBox=0;

figure
Z=lscan;
Z(Z<=0)=.05;
Z=log10(Z);
[a,b]=size(lscan);

[l1,l2]=size(lscan);
xb=22.15;%-.008333333;%;%xbottom
yb=112.789;%-.008333333;%+30*.008333333;

x=locs(1:end-3,2);%2:
x=ceil((x-xb)/.008333333);
x=a+1-x;%165-x;
y=locs(1:end-3,3);%2:
y=ceil((y-yb)/.008333333);
f=[min(x),max(x),min(y),max(y)];
for i=1:length(x)
    Z(x(i),y(i))=-2;
end
%}
minx=min(x); maxx=max(x); miny=min(y); maxy=max(y);
[l1,l2]=size(Z);
%Z(Z<0)=NaN;
Z=Z(minx:maxx+4,miny-4:maxy+4);

Z2=lscan(minx:maxx+4,miny-4:maxy+4);
y=y-miny+5;

if truncateBox==1
    [z1,z2]=size(Z);
    lx=41;%x-length
    yto=63;%y-length
    xfrom=z1-lx+1;
    xdiff=z1-lx;
    findx=find(x<xfrom);
    x(findx)=[]; y(findx)=[];
    findy=find(y>yto);
    x(findy)=[]; y(findy)=[];
    x=x-xdiff;
    Z=Z(xfrom:end,1:yto);
    Z2=Z2(xfrom:end,1:yto);
end

colormap([gray;.8,.0,.0])
%colormap([hot;0,.25,.5])
colormap(flipud(colormap))
la=50; lb=150;
imagesc(Z)
%axis([miny-.5,maxy+.5,minx-.5,maxx+.5])
xlabel('Latitude')
ylabel('Longitude')
set(gca,'xtick',[],'ytick',[])