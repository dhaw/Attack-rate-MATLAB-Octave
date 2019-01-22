function f=plotinspace(D,fglob)%,minNind,maxNind)
%trim=4; trim2=2*trim;
%%D=reshape(D,33-trim2,55-trim2);%Or trimbyk then reshape
%%[a,b]=size(D);

%D=D(trim+1:end-trim,trim-1:end-trim);

maxD=max(max(D));
figure
fs=12;
%colormap parula
%colormap gray
%numCol=20;
%cmap=flipud(gray(numCol));
%cmap=colormap(gray);
%cmap=flipud(cmap);

%cmap=flipud(gray(100));
cmap=redblue(100);
D=D-fglob;

%cmap=[1,1,1;.5,.5,.5];
%cmap=cmap.*repmat([.165,.31,.431]/.431,numCol,1);
%cmap=cmap(numCol/2+1:end,:);
%cmap=othercolor('Oranges9');
%cmap=parula(50);
%cmap(1,:)=[1,1,1];
%cmap=[1,1,1;cmap];
colormap(cmap);
%D(D==-inf)=0;

%imx=zeros(a+trim2,b+trim2);
%imx(minNind)=1; imx(maxNind)=1;
%imx=imx(trim+1:end-trim,trim+1:end-trim);
%imx=flipud(imx);
%vtop=maxD;%max(1,max(max(D)));%maxD;%floor(maxD/.05);

vtop=max(max(D));
vbot=min(min(D));
vmax=max(abs([vbot,vtop]));

vv=floor(100*vmax)/10;

imagesc(D);
%im2=imagesc(imx);
%set(im2,'AlphaData',.5)

xlabel('Longitude','fontsize',fs)
ylabel('Latitude','fontsize',fs)
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
set(gca,'fontsize',fs)
hcb=colorbar;

caxis([-vmax,vmax])%([0,vmax])([-vmax,vmax])%cmap=[1,1,1;cmap];

v=(-vv:.02:vv);
%set(hcb,'YTick',v)
set(hcb,'YTick',-.4:.1:.4)
hold off

f=maxD;
%{
D=flipud(D);
[a,b]=size(D);
maxD=max(max(D));
figure
colormap parula
%colormap gray
%cmap=gray(100);
%cmap(1,:)=[1,1,1];
%cmap=[1,1,1;cmap];
%colormap(cmap);
%D(D==-inf)=-1;

imx=ones(a,b);
imx(minNind)=0; imx(maxNind)=0;
imx=flipud(imx);

figure
hold on
im1=imagesc(D);
im2=imagesc(imx);
set(im2,'AlphaData',.5)
xlabel('longitude')
ylabel('latitude')
hcb=colorbar;
vtop=1;%maxD;%floor(maxD/.05);
caxis([0 vtop])%cmap=[1,1,1;cmap];

v=(0:.05:.05*vtop);
%set(hcb,'YTick',v)
f=maxD;
%}