function [f,g]=aggResult(D,AR,k)%,NN
[a,b]=size(D);
A=reshape(D,a*b,1);
NN=A;
meanA1=nansum(AR(:,1).*NN)./sum(NN);
Anan=1-isnan(A);
A(Anan==0)=0;
meanA=sum(A.*NN)./sum(Anan.*NN);
n=a*b;
%x=nanmean(AR(:,6:end),2);
x=AR;
x=x.*NN;
x=reshape(x,a,b);
N=D;%reshape(D,a,b);
ak=k*floor(a/k);
bk=k*floor(b/k);
x=x(a-ak+1:end,1:bk);
N=N(a-ak+1:end,1:bk);
D=N; a=ak; b=bk; n=ak*bk;
%{
lH=length(IH);
X2=zeros(n/k^2,lH);%NN2=X2;
%Y=X;
x2=zeros(a/k,b/k);n2=x2;
for i=1:lH
    H=cell2mat(IH(i));
    x=sum(H>0,2);
    x=reshape(x,a,b);
    %y=sum(H>0,2);
    %
%}
    for k1=1:k:a-k+1
        idxr=1+floor(k1/k);
        for k2=1:k:b-k+1
            idxc=1+floor(k2/k);
            MF=x(k1:k1+k-1,k2:k2+k-1);
            NF=D(k1:k1+k-1,k2:k2+k-1);
            tempn=max(sum(NF(:)),1);% tempn(tempn==0)=1;
            n2(idxr,idxc)=tempn;
            x2(idxr,idxc)=sum(MF(:))/n2(idxr,idxc);
        end
    end
    X2=reshape(x2,n/k^2,1);
    NN2=reshape(n2,n/k^2,1); NN20=NN2; NN20(NN2==0)=1;
%end

f=X2;%var(X2);
g=NN2;
%
logN2=log10(NN2);
maxN=max(NN2);%max(logN2);
%{
figure
fs=40;

semilogx(NN2,X2,'x','color',[.447,.553,.647],'markersize',30,'LineWidth',3);%[.165,.31,.431][.447,.553,.647]
hold on
semilogx([1,maxN],[meanA1,meanA1],'k--','linewidth',3);

%plot([0,maxN],[meanA,meanA],'k--','linewidth',1);
minY=min(min(X2));
maxY=max(max(X2));
%plot(logN2,X2,'o','color',[0,.5,.5]);
xlabel('log N','FontSize',fs);
ylabel('Proportion immune','FontSize',fs);
set(gca,'FontSize',fs);
minY=max(minY-.1,0); maxY=min(maxY+.1,1);
axis ([0,maxN,0,maxY])%minY,maxY])
%legend('Rural','Urban')
grid on
hold off
%}
end
    %{
figure
M=[X2,NN2];
M=sortrows(M,lH+1);
M(:,end)=[];
M=flipud(M);
cmap=othercolor('Oranges9');
colormap(cmap)
imagesc(1:lH,1:n/k^2,M)
%set(gca,'XTick',(10:10:lH),'YTick',(5:5:n))
hcb=colorbar;
caxis([0 1])
set(hcb,'YTick',(0:.2:1))
xlabel('Time \tau')
ylabel('Cell')
title('Proportion of immune agents')
    %}
%{
%15x15:
g=2;
%
%10x10:
g=3;
%
%5x5:
g=6;
%
%}