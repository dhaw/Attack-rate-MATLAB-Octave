function [f,g]=selectSubgrid(D,x,y)
[a,b]=size(D);
X=zeros(a,b);
xy=sub2ind([a,b],x,y);
X(xy)=1;
ncells=length(x);
grid=X;
%minx=zeros(a,1); maxx=minx;
%miny=zeros(b,1); maxy=miny;
%{
for i=1:a
    Xi=X(i,:);
    findx=find(Xi==1);
    if isempty(findx)==0
        minx=max(findx(1)-4,1);
        maxx=min(findx(end)+4,a);
        grid(i,minx:maxx)=1;
    end
end
for j=1:b
    Xj=X(:,j);
    findy=find(Xj==1);
    if isempty(findy)==0
        miny=max(findy(1)-4,1);
        maxy=min(findy(end)+4,b);
        grid(miny:maxy,j)=1;
    end
end
%}
xm=max(x-4,1);
xp=min(x+4,a);
ym=max(y-4,1);
yp=min(y+4,b);
for i=1:ncells
    xmi=xm(i); xpi=xp(i);
    ymi=ym(i); ypi=yp(i);
    grid(xmi:xpi,ymi:ypi)=1;
end
    
grid2=grid;

for i=1:a
    Xi=grid(i,:);
    findx=find(Xi==1);
    if isempty(findx)==0
        minx=findx(1);
        maxx=findx(end);
        grid2(i,minx:maxx)=1;
    end
end
for j=1:b
    Xj=grid(:,j);
    findy=find(Xj==1);
    if isempty(findy)==0
        miny=findy(1);
        maxy=findy(end);
        grid2(miny:maxy,j)=1;
    end
end
%{
figure
grid3=grid2;
grid3(xy)=2;
imagesc(grid3);
%}
%
grid4=-grid2;
grid4(xy)=1:ncells;
XY=reshape(grid4,a*b,1);
XY(XY==0)=[];
XY(XY<0)=0;
f=grid2;%1 in cells to keep
g=XY;%fluscape location numbers
%}


%Old first step - upper y border didn't work, also had cells too close to
%border
%f=[a*b,sum(sum(grid2))];
%{
xfx=xf;
yfx=yf;
xfx(xf<10)=NaN; yfx(xf<10)=NaN;
xfx(yf>65)=NaN; yfx(yf>65)=NaN;
xfx(isnan(xfx)==1)=[]; yfx(isnan(yfx)==1)=[];
selectSubgrid(Df,xfx,yfx)
%}
