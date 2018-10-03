function [Z2,x,y]=transectplot(lscan,transect,locs)%Df,Dt,locs
%figure
Z=lscan;
%Z(Z<=0)=.05;
Z=log10(Z);
ZZ=Z;
Z(isnan(Z)==1)=0;
Z(isinf(Z)==1)=0;
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
ZZ=ZZ(minx:maxx+4,miny-4:maxy+4);
minxaxis=min(locs(:,2);

Z2=lscan(minx:maxx+4,miny-4:maxy+4);
y=y-miny+5;
%{
%Z(transect==0)=NaN;
%transect2=transect;
%transect2(1:9,:)=0; transect2(:,66:end)=0;
%bigMask=imdilate(transect, true(1));
colormap([gray;.8,.0,.0])
%colormap([hot;0,.25,.5])
colormap(flipud(colormap))
la=50; lb=150;
imagesc(Z)
hold on
%visboundaries(transect2,'linestyle','--','color',[0,.8,.8])
visboundaries(transect,'linestyle','-','linewidth',2,'EnhanceVisibility',false,'color',[0,.8,.8])
colorbar
%}

Z(Z~=-2)=NaN; Z(Z==-2)=1;
Y=transect;
Y(Y~=1)=0; Y=1-Y;
%Y=flipud(1-Y);
fs=12;
figure
ax1=axes;
imagesc(ax1,ZZ);
xlabel('Longitude')
ylabel('Latitude')
set(gca,'xtick',[],'ytick',[],'fontsize',fs)
ax2=axes;
b=imagesc(ax2,Z);
set(b,'AlphaData',~isnan(Z));
set(gca,'xtick',[],'ytick',[],'fontsize',fs)
ax3=axes;
%c=visboundaries(Y,'linestyle','-','linewidth',2,'EnhanceVisibility',false,'color',[0,.8,.8]);%
c=imagesc(ax3,Y,'AlphaData',.3*Y);
%set(c,'AlphaData',isnan(Y))
set(gca,'xtick',[],'ytick',[],'fontsize',fs)
linkaxes([ax1,ax2,ax3])
set([ax1,ax2,ax3],'Position',[.1 .11 .685 .815]);%[.17 .11 .685 .815]); ,ax3
ax2.Visible='off';
ax3.Visible='off';
cmap=colormap(gray);
colormap(ax1,flipud(cmap))
colormap(ax2,[1,0,0])
colormap(ax3,[.5,.5,.5,;0,0,1])%[0,1,1,;1,1,1]
colorbar(ax1,'Position',[.82 .11 .05 .815],'fontsize',fs);%[.88 .11 .0675 .815],'fontsize',fs);
%{
colormap([gray])
cmap=flipud(colormap);
colormap(cmap);
a=imagesc(ZZ);
%colorbar
hold on
Z(Z~=-2)=NaN; Z(Z==-2)=1;
%colormap([0.8,0,0])
b=imagesc(Z);
set(b,'AlphaData',~isnan(Z))
Y=transect;
Y(Y~=1)=NaN;
%colormap([0,0,0;0,.8,.8])
c=imagesc(Y,'alphadata',.1);
%set(c,'AlphaData',.5);
set(c,'AlphaData',isnan(Y))
%}

%axis([miny-.5,maxy+.5,minx-.5,maxx+.5])
xlabel('Latitude')
ylabel('Longitude')
set(gca,'xtick',[],'ytick',[],'fontsize',fs)