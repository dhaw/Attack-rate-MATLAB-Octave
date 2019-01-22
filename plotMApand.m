function f=plotMApand(A,Acum,NN,fp)
%A and Acum/cumulative - wring way round! A=f, Acum=g
cbaroff=1;
fs=12;%Font size: 30 for large
lw=2;%1.5
tauend=size(A,2);
newend=min(tauend,500);
cell=6;
Ni=repmat(NN,1,tauend);
B=sum(Acum.*Ni,1)/sum(NN);
n=length(NN);
%%
%As fraction of total AR for each year:
%{
sumYear=sum(A,1);
sumYear(sumYear==0)=1;
sumYear=repmat(sumYear,n,1);
A=A./sumYear;
sumYear=sum(Acum,1);
sumYear(sumYear==0)=1;
sumYear=repmat(sumYear,n,1);
Acum=Acum./sumYear;
%}
%%
%Trim:
%{
NN=trimbyk(NN); n=length(NN); A2=zeros(n,tauend); A3=A2;
for i=1:tauend
    A2(:,i)=trimbyk(Acum(:,i));
    A3(:,i)=trimbyk(A(:,i));
end
Acum=A2; A=A3;
%}
%%
%cmap=[cool(n);gray(n)];
NA=[NN,Acum]; NA=sortrows(NA,1); Acum=NA(:,2:end);
NA=[NN,A]; NA=sortrows(NA,1); A=NA(:,2:end);
Nsort=sortrows(NN); [maxN,cellm]=max(Nsort);
%cellm=200;
%f=Nsort(cellm);
%%
figure
%hold on
logN=log10(Nsort); logN(Nsort==0)=0;
cmap=colormap(parula);
%cmap=.9*cmap;
colormap(cmap)
cc=colormap;
%imagesc(logN);
%cc=flipud(colormap);
lc=size(cc,1);
if n>1
    Gs=round(interp1(linspace(min(logN(:)),max(logN(:)),lc),1:lc,logN));
    CC=reshape(cc(Gs,:),[size(Gs) 3]);%Make RGB image from scaled
    %
    %thismany=20; sam=randsample(1:n,thismany);
    %sam=(600+rem(n,50):50:n); thismany=length(sam);
    nn=1;
    bottom=rem(n,nn);
    if bottom==0
        bottom=1;
    end
    sam=(bottom:nn:n); thismany=length(sam);%nn+
    y1=Acum(sam,:); y2=A(sam,:); cc=CC(sam,:);
else
    y1=Acum; y2=A; thismany=1; cc=[.447,.553,.647];
end
%
y1=y1(:,1:newend); y2=y2(:,1:newend);
T=1:newend;

X=zeros(1,newend);
Y=X;
for i=1:newend
    xi=corrcoef(Acum(:,i),fp);
    X(i)=xi(1,2);
    yi=corrcoef(A(:,i),fp);
    Y(i)=yi(1,2);
end

hold on
for i=1:thismany%Yes, this loop is clonky
    plot(T,y1(i,:),':','linewidth',lw,'color',cc(i,:));%[.5,.5,.5])%cc(i,:)); 'o-'
    plot(T,y2(i,:),'-','linewidth',lw,'color',cc(i,:));%[0,0,0])%cc(i,:)); 'o-'
end
plot(T,X,'k-','linewidth',lw)
plot(T,Y,'k--','linewidth',lw)
xlabel('Time (years)','FontSize',fs)
ylabel('Attack rate')%('Proportion immune','FontSize',fs) Relative attack
set(gca,'FontSize',fs);
maxY=max(max([y1;y2])); maxY=min(maxY+.1,1);%+.1
axis tight%([0,newend,0,maxY])
grid on
grid minor
box on
if n>1 && cbaroff~=1
    colorbar
    caxis([0,max(logN)])
end
hold off
%}
%{
figure
%colormap=cmap;
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