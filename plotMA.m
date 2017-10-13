function f=plotMA(A,Acum,NN)
fs=15;%Font size: 30 for large
tauend=size(A,2);
cell=6;
Ni=repmat(NN,1,tauend);
B=sum(Acum.*Ni,1)/sum(NN);
n=length(NN);
%Trim:
%{
NN=trimbyk(NN); n=length(NN); A2=zeros(n,tauend); A3=A2;
for i=1:tauend
    A2(:,i)=trimbyk(Acum(:,i));
    A3(:,i)=trimbyk(A(:,i));
end
Acum=A2; A=A3;
%}

%cmap=[cool(n);gray(n)];
NA=[NN,Acum]; NA=sortrows(NA,1); Acum=NA(:,2:end);
NA=[NN,A]; NA=sortrows(NA,1); A=NA(:,2:end);
Nsort=sortrows(NN); [maxN,cellm]=max(Nsort);
%cellm=200;
%f=Nsort(cellm);

%
figure
%hold on
logN=log10(Nsort); logN(Nsort==0)=0;
%cmap=othercolor('BrBG8');
cmap=colormap(parula);
%cmap=.9*cmap;
colormap(cmap)
%imagesc(logN);
cc=colormap;
%cc=flipud(colormap);
lc=size(cc,1);
if n>1
    Gs=round(interp1(linspace(min(logN(:)),max(logN(:)),lc),1:lc,logN));
    CC=reshape(cc(Gs,:),[size(Gs) 3]);%Make RGB image from scaled
    %hold off
    %
    %thismany=20; sam=randsample(1:n,thismany);
    nn=2;
    %sam=(600+rem(n,50):50:n); thismany=length(sam);
    sam=(nn+rem(n,nn):nn:n); thismany=length(sam);
    y1=Acum(sam,:); y2=A(sam,:); cc=CC(sam,:);
    %figure
else
    y1=Acum; y2=A; thismany=1; cc=[.447,.553,.647];
end
T=1:tauend;
hold on
for i=1:thismany%Yes, this loop is clonky
    plot(T,y1(i,:),'o-','linewidth',1.5,'color',cc(i,:));
    plot(T,y2(i,:),'o:','linewidth',1.5,'color',cc(i,:));
end
xlabel('Time (years)','FontSize',fs)
ylabel('Proportion immune','FontSize',fs)
set(gca,'FontSize',fs);
maxY=max(max([y1;y2])); maxY=min(maxY+.1,1);%+.1
axis([1,tauend,0,maxY])
grid on
if n>1
    colorbar
    caxis([0,max(logN)])
end
hold off
%}
%{
figure
colormap=cmap;
hold on
h1=plot(1:tauend,Acum,'-o','linewidth',2,'markerfacecolor','auto');
h2=plot(1:tauend,A,':o','linewidth',2,'markerfacecolor','auto');
%h3=plot(1:tauend,Acum(cellm,:),'k-o','linewidth',2,'markerfacecolor','auto');
%h4=plot(1:tauend,A(cellm,:),'k:o','linewidth',2,'markerfacecolor','auto');
%plot(1:tauend,B,'k-o','linewidth',1.5,'markerfacecolor','k');
xlabel('Time','FontSize',fs)
ylabel('Proportion immune','FontSize',fs)
set(gca,'FontSize',fs);
axis([1,tauend,0,1])
grid on
hold off
%}