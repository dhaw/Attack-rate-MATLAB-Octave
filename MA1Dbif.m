function [f,g]=MA1Dbif%(X,X2)
eps=.13;
cross=1;
param=(0:.01:1);%CROSS in MA1D2sub
%
tauend=1000;%Must match MA1D2sub
lp=length(param);
burn=200;
B=-1*ones(lp,tauend-burn);
B2=B;
numPoints=zeros(lp,1);
numPoints2=numPoints;
for i=1:lp
    parami=param(i);
    %{
    %Cross:
    [ignore,g]=MA1D2sub(eps,parami);
    g=g(:,burn+1:end);
    gsum=sum(g,1); g1=g(1,:);
    G=g1./gsum; G(gsum==0)=0;
    Gunique=unique(G);
    lGu=length(Gunique);
    B(i,1:lGu)=Gunique;
    numPoints(i)=lGu;
    %}
    %
    %Eps:
    %[f,g]=MA1D2sub(parami,cross);
    [f,g,other]=MA1D(2,0,parami);
    f=f(:,burn+1:end);
    %f=1-f(:,1);%Out if use MA1D
    Funique=unique(f);
    lFu=length(Funique);
    B2(i,1:lFu)=Funique;
    numPoints2(i)=lFu;g;%
    %
    g=g(:,burn+1:end);
    G=g;%sum(g,1);%G=g if use MA1D
    Gunique=unique(G);
    lGu=length(Gunique);
    B(i,1:lGu)=Gunique;
    numPoints(i)=lGu;
    %}
end
%
numPoints=max(numPoints);
X=B(:,1:numPoints);
f=X;
%
%Eps only:
numPoints2=max(numPoints2);
X2=B2(:,1:numPoints2);
%
g=X2;
%}
%
fs=12; lw=2; ms=2.5;%7.5;
figure
%{
%Cross:
plot(param,X','k.','markersize',ms)
axis([param(1),param(end),0,1])
xlabel('Cross immunity \chi')
ylabel('R^1/(R^1+R^2)')
set(gca,'fontsize',fs)
grid on
grid minor
box on
%}
%
%Eps:
hold on
plot(param,X2','.','markersize',ms,'color',[0,0,0])
plot(param,X','.','markersize',ms,'color',[0,0,.8])
hold off
axis([param(1),param(end),0,.4])
xlabel('\epsilon')
ylabel('Proportion immune')
set(gca,'fontsize',fs)
%legend('Total immune','Newly infected')
grid on
grid minor
box on
%}

    