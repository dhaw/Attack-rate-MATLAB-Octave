function f=plotARDist(f1,f2,a,b)
%{
Btot=B1+B2+B3+B4; Bmean=mean(mean(Btot));
Bsum=sum(sum(Btot));
B1p=B1./Btot; B1p=B1p-sum(sum(B1))/Bsum;%B1p=B1p-sum(sum(B1))/sum(sum(Btot));
B2p=B2./Btot; B2p=B2p-sum(sum(B2))/Bsum;%B2=B2-mean(mean(B2));
B3p=B3./Btot; B3p=B3p-sum(sum(B3))/Bsum;%B3=B3-mean(mean(B3));
B4p=B4./Btot; B4p=B4p-sum(sum(B4))/Bsum;%B4p=B4p-mean(mean(B4p));
%}
A1=reshape(f1,a,b);
A2=reshape(f2,a,b);
minPop=0;
maxPop=(max([f1;f2]));
%limit=max(abs([minPop,maxPop]));
%minPop=-limit; maxPop=limit;

figure

%Using panel:
p=panel();
p.pack({1}, {.45,.45,.1});%(1,2);
%p(1,1).pack({.5,.5});

% set margins
p.de.margin = 10;
%%p(1,1).marginbottom = 20;
%p(1,2).marginleft = 0;
p.margin = [5,5,5,5];

% set some font properties
p.select('all');
p.fontname = 'Helvetica';
p.fontsize = 12;
p.fontweight = 'bold';

%p.identify();

cmap=parula(100);%redblue
colormap(cmap)

p(1,1).select();
imagesc(A1)
axis image
caxis([minPop,maxPop])
%cbar=colorbar;
title('Age-independent kernel')
set(gca,'xtick',[])
set(gca,'ytick',[])

p(1,2).select();
imagesc(A2)
axis image
caxis([minPop,maxPop])
title('Age-dependent kernel')
set(gca,'xtick',[])
set(gca,'ytick',[])

p(1,3).select();
c=colorbar;
x=get(c,'Position');
x(3)=0.05;
x(1)=.89;
set(c,'Position',x)
axis off
caxis([minPop,maxPop])
set(gca,'xtick',[])
set(gca,'ytick',[])