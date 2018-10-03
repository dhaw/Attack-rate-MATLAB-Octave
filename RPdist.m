function f=RPdist(tout,yout,boxLat,boxLong)
%Grid size:
a=boxLat;
b=boxLong;

[lt,ncells]=size(yout);
peakInd=zeros(1,ncells);
peakTime=peakInd;
for i=1:ncells
    yi=yout(:,i);
    [maxi,findi]=max(yi);
    findi=findi(1);
    peakInd(i)=findi;
    peakTime(i)=tout(findi);
end
X=reshape(peakTime,a,b);
figure
fs=18;
colormap parula
imagesc(X)
set(gca,'fontsize',fs,'xticklabels',[],'yticklabels','')



