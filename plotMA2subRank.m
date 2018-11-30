function f=plotMA2subRank(g)%,NN)
[a,tauend]=size(g);
nbar=a/2;
n=nbar/4;
h1=g(1:n,:)+g(n+1:2*n,:)+g(2*n+1:3*n,:)+g(3*n+1:4*n,:);
h2=g(4*n+1:5*n,:)+g(5*n+1:6*n,:)+g(6*n+1:7*n,:)+g(7*n+1:end,:);
h1sum=sum(h1,1); h2sum=sum(h2,1);
h1=h1(:,h1sum>.0001); h2=h2(:,h2sum>.0001);
l1=size(h1,2); l2=size(h2,2);
X1=zeros(n,l1); X2=zeros(n,l2);
nvec=(1:n)';
for i=1:l1
    v=[nvec,h1(:,i)];
    v=sortrows(v,2);
    X1(:,i)=v(:,1);
end
for i=1:l2
    v=[nvec,h2(:,i)];
    v=sortrows(v,2);
    X2(:,i)=v(:,1);
end
fs=15; lw=2;
figure
    subplot(2,1,1)
    plot(1:l1,X1,'-','linewidth',lw)%,'color',[.5,0,0])%'o-'
    xlabel('Time (years)','FontSize',fs)
    ylabel('Rank (H1N1)','FontSize',fs)
    axis ([0,l1,0,n]);
    axis tight
    set(gca,'FontSize',fs);
    grid on
    grid minor
    box on
    subplot(2,1,2)
    plot(1:l2,X2,'-','linewidth',lw)%,'color',[0,0,.5])%'o-'
    xlabel('Time (years)','FontSize',fs)
    ylabel('Rank (H3N2)','FontSize',fs)
    axis ([0,l2,0,n]);
    axis tight
    set(gca,'FontSize',fs);
    grid on
    grid minor
    box on

%{
asum1=sum(F,1);
asum2=sum(G,1);
asum=asum1.*asum2;
Arank=Arank(:,asum>.001);
AcumRank=AcumRank(:,asum>.001);

cbaroff=1;
fs=15;%Font size: 30 for large
lw=2;%1.5
tauend=size(Arank,2);
newend=min(tauend,500);
Ni=repmat(NN,1,tauend);
n=length(NN);

figure
%hold on
colormap(parula);
T=1:newend;
plot(T,Arank,'-')%,'linewidth',lw,'color',cc(i,:));%[0,0,0])%cc(i,:)); 'o-'
xlabel('Time (years)','FontSize',fs)
ylabel('Rank')%('Proportion immune','FontSize',fs) Relative attack
set(gca,'FontSize',fs);
%maxY=max(max([y1;y2])); maxY=min(maxY+.1,1);%+.1
axis([0,newend,0,n])
grid on
grid minor
box on

figure
plot(T,AcumRank,'-')%,'linewidth',lw,'color',cc(i,:));%[0,0,0])%cc(i,:)); 'o-'
xlabel('Time (years)','FontSize',fs)
ylabel('Attack rate')%('Proportion immune','FontSize',fs) Relative attack
set(gca,'FontSize',fs);
%maxY=max(max([y1;y2])); maxY=min(maxY+.1,1);%+.1
axis([0,newend,0,n])
grid on
grid minor
box on
%}
%imagesc(logN);
%cc=flipud(colormap);
%{
lc=size(cc,1);
if n>1
    Gs=round(interp1(linspace(min(logN(:)),max(logN(:)),lc),1:lc,logN));
    CC=reshape(cc(Gs,:),[size(Gs) 3]);%Make RGB image from scaled
    %
    %thismany=20; sam=randsample(1:n,thismany);
    %sam=(600+rem(n,50):50:n); thismany=length(sam);
    nn=10;
    bottom=rem(n,nn);
    if bottom==0
        bottom=1;
    end
    sam=(bottom:nn:n); thismany=length(sam);%nn+
    y1=Ar(sam,:); y2=A(sam,:); cc=CC(sam,:);
else
    y1=Acum; y2=A; thismany=1; cc=[.447,.553,.647];
end
%
y1=y1(:,1:newend); y2=y2(:,1:newend);
T=1:newend;
hold on
for i=1:thismany%Yes, this loop is clonky
    plot(T,y1(i,:),':','linewidth',lw,'color',cc(i,:));%[.5,.5,.5])%cc(i,:)); 'o-'
    plot(T,y2(i,:),'-','linewidth',lw,'color',cc(i,:));%[0,0,0])%cc(i,:)); 'o-'
end
xlabel('Time (years)','FontSize',fs)
ylabel('Attack rate')%('Proportion immune','FontSize',fs) Relative attack
set(gca,'FontSize',fs);
maxY=max(max([y1;y2])); maxY=min(maxY+.1,1);%+.1
axis([0,newend,0,maxY])
grid on
grid minor
if n>1 && cbaroff~=1
    colorbar
    caxis([0,max(logN)])
end
hold off
%}
