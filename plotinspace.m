function f=plotinspace(D,minNind,maxNind)
trim=4; trim2=2*trim;
D=reshape(D,33-trim2,55-trim2);%Or trimbyk then reshape
[a,b]=size(D);

maxD=max(max(D));
figure
fs=15;
%colormap parula
%colormap gray
%cmap=flipud(gray(20));
%cmap=othercolor('Oranges9');
cmap=parula(50);
%cmap(1,:)=[1,1,1];
cmap=[cmap;0,0,0];
colormap(cmap);
%D(D==-inf)=0;

imx=zeros(a+trim2,b+trim2);
imx(minNind)=1; imx(maxNind)=1;
imx=imx(trim+1:end-trim,trim+1:end-trim);
imx=flipud(imx);
vtop=maxD;%max(1,max(max(D)));%maxD;%floor(maxD/.05);

imagesc(D);
%im2=imagesc(imx);
%set(im2,'AlphaData',.5)
xlabel('longitude','fontsize',fs)
ylabel('latitude','fontsize',fs)
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
set(gca,'fontsize',fs)
hcb=colorbar;
caxis([0 vtop])%cmap=[1,1,1;cmap];
v=(0:.05:.05*vtop);
set(hcb,'YTick',v)
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