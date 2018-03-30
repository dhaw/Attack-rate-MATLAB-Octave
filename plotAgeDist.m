<<<<<<< HEAD
function f=plotAgeDist(B1,B2,B3,B4)
isPop=1;
%
Btot=B1+B2+B3+B4; Bmean=mean(mean(Btot));
Bsum=sum(sum(Btot));
B1p=B1./Btot; B1p=B1p-sum(sum(B1))/Bsum;%B1p=B1p-sum(sum(B1))/sum(sum(Btot));
B2p=B2./Btot; B2p=B2p-sum(sum(B2))/Bsum;%B2=B2-mean(mean(B2));
B3p=B3./Btot; B3p=B3p-sum(sum(B3))/Bsum;%B3=B3-mean(mean(B3));
B4p=B4./Btot; B4p=B4p-sum(sum(B4))/Bsum;%B4p=B4p-mean(mean(B4p));
%}
minPop=(min(min([B1p;B2p;B3p;B4p])));
maxPop=(max(max([B1p;B2p;B3p;B4p])));%log10
limit=max(abs([minPop,maxPop]));
minPop=-limit; maxPop=limit;

%Using panel:
p=panel();
p.pack({1}, {.9,.1});%(1,2);
p(1,1).pack({.5,.5}, {.5,.5});

% set margins
p.de.margin = 10;
%p(1,1).marginbottom = 20;
p(1,2).marginleft = 0;
p.margin = [5,5,5,5];

% set some font properties
p.select('all');
p.fontname = 'Helvetica';
p.fontsize = 12;
p.fontweight = 'bold';

%p.identify();

cmap=parula(100);%redblue
colormap(cmap)

p(1,1,1,1).select();
imagesc(B1p)
axis image
caxis([minPop,maxPop])
%cbar=colorbar;
title('0-4')
set(gca,'xtick',[])
set(gca,'ytick',[])

p(1,1,1,2).select();
imagesc(B2p)
axis image
caxis([minPop,maxPop])
title('5-19')
set(gca,'xtick',[])
set(gca,'ytick',[])

p(1,1,2,1).select();
imagesc(B3p)
axis image
caxis([minPop,maxPop])
title('20-64')
set(gca,'xtick',[])
set(gca,'ytick',[])

p(1,1,2,2).select();
imagesc(B4p)
axis image
caxis([minPop,maxPop])
title('64+')
set(gca,'xtick',[])
set(gca,'ytick',[])

p(1,2).select();
c=colorbar;
x=get(c,'Position');
x(3)=0.05;
x(1)=.89;
set(c,'Position',x)
axis off
caxis([minPop,maxPop])
set(gca,'xtick',[])
set(gca,'ytick',[])

%{
figure
%suptitle('Age distribution')
cmap=redblue(100);
colormap(cmap)
subplot(2,2,1)
imagesc(B1p)
axis image
caxis([minPop,maxPop])
title('0-4')
set(gca,'xtick',[])
set(gca,'ytick',[])
%
subplot(2,2,2)
imagesc(B2p)
axis image
caxis([minPop,maxPop])
title('5-19')
set(gca,'xtick',[])
set(gca,'ytick',[])
%
subplot(2,2,3)
imagesc(B3p)
axis image
caxis([minPop,maxPop])
title('20-64')
set(gca,'xtick',[])
set(gca,'ytick',[])
%
subplot(2,2,4)
imagesc(B4p)
axis image
caxis([minPop,maxPop])
title('65+')
set(gca,'xtick',[])
set(gca,'ytick',[])

hp4=get(subplot(2,2,4),'Position');
colorbar('Position',[hp4(1)+hp4(3)+0.01 ,hp4(2),0.02,hp4(2)+hp4(3)*2.1])
=======
function f=plotAgeDist(B1,B2,B3,B4)
%
Btot=B1+B2+B3+B4; Bmean=mean(mean(Btot));
Bsum=sum(sum(Btot));
B1p=B1./Btot; B1p=B1p-sum(sum(B1))/Bsum;%B1p=B1p-sum(sum(B1))/sum(sum(Btot));
B2p=B2./Btot; B2p=B2p-sum(sum(B2))/Bsum;%B2=B2-mean(mean(B2));
B3p=B3./Btot; B3p=B3p-sum(sum(B3))/Bsum;%B3=B3-mean(mean(B3));
B4p=B4./Btot; B4p=B4p-sum(sum(B4))/Bsum;%B4p=B4p-mean(mean(B4p));
%}
minPop=(min(min([B1p;B2p;B3p;B4p])));
maxPop=(max(max([B1p;B2p;B3p;B4p])));%log10
limit=max(abs([minPop,maxPop]));
minPop=-limit; maxPop=limit;

%Using panel:
p=panel();
p.pack({1}, {.9,.1});%(1,2);
p(1,1).pack({.5,.5}, {.5,.5});

% set margins
p.de.margin = 10;
%p(1,1).marginbottom = 20;
p(1,2).marginleft = 0;
p.margin = [5,5,5,5];

% set some font properties
p.select('all');
p.fontname = 'Helvetica';
p.fontsize = 12;
p.fontweight = 'bold';

%p.identify();

cmap=redblue(100);
colormap(cmap)

p(1,1,1,1).select();
imagesc(B1p)
axis image
caxis([minPop,maxPop])
%cbar=colorbar;
title('0-4')
set(gca,'xtick',[])
set(gca,'ytick',[])

p(1,1,1,2).select();
imagesc(B2p)
axis image
caxis([minPop,maxPop])
title('5-19')
set(gca,'xtick',[])
set(gca,'ytick',[])

p(1,1,2,1).select();
imagesc(B3p)
axis image
caxis([minPop,maxPop])
title('20-64')
set(gca,'xtick',[])
set(gca,'ytick',[])

p(1,1,2,2).select();
imagesc(B4p)
axis image
caxis([minPop,maxPop])
title('64+')
set(gca,'xtick',[])
set(gca,'ytick',[])

p(1,2).select();
c=colorbar;
x=get(c,'Position');
x(3)=0.05;
x(1)=.89;
set(c,'Position',x)
axis off
caxis([minPop,maxPop])
set(gca,'xtick',[])
set(gca,'ytick',[])

%{
figure
%suptitle('Age distribution')
cmap=redblue(100);
colormap(cmap)
subplot(2,2,1)
imagesc(B1p)
axis image
caxis([minPop,maxPop])
title('0-4')
set(gca,'xtick',[])
set(gca,'ytick',[])
%
subplot(2,2,2)
imagesc(B2p)
axis image
caxis([minPop,maxPop])
title('5-19')
set(gca,'xtick',[])
set(gca,'ytick',[])
%
subplot(2,2,3)
imagesc(B3p)
axis image
caxis([minPop,maxPop])
title('20-64')
set(gca,'xtick',[])
set(gca,'ytick',[])
%
subplot(2,2,4)
imagesc(B4p)
axis image
caxis([minPop,maxPop])
title('65+')
set(gca,'xtick',[])
set(gca,'ytick',[])

hp4=get(subplot(2,2,4),'Position');
colorbar('Position',[hp4(1)+hp4(3)+0.01 ,hp4(2),0.02,hp4(2)+hp4(3)*2.1])
>>>>>>> acfad3e4da850e525286e08138beef9f06a8c2cc
%}