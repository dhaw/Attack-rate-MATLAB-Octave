function f=MA1Dbif%(X)
eps=.1;
cross=(0:.01:1);%CROSS in MA1D2sub
%
tauend=500;%Must match MA1D2sub
lc=length(cross);
burn=200;
B=-1*ones(lc,tauend-burn);
numPoints=zeros(lc,1);
for i=1:lc
    ci=cross(i);
    [ignore,g]=MA1D2sub(eps,ci);
    g=g(:,burn+1:end);
    gsum=sum(g,1); g1=g(1,:);
    G=g1./gsum; G(gsum==0)=0;
    Gunique=unique(G);
    lGu=length(Gunique);
    B(i,1:lGu)=Gunique;
    numPoints(i)=lGu;
end
numPoints=max(numPoints);
X=B(:,1:numPoints);
%}
f=X;
%
fs=15; lw=2; ms=2;
figure
plot(cross,X','ko','markersize',ms)
axis([cross(1),cross(end),0,1])
xlabel('Cross immunity')
ylabel('R^1/(R^1+R^2)')
set(gca,'fontsize',fs)
grid on
grid minor

    