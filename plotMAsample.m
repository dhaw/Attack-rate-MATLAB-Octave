function f=plotMAsample(A,Acum,NN,sample)
%sample=find(Is)
%For now: colormap not rep of pop density - relevant for fluscape?
cbaroff=1;
fs=12;%Font size: 30 for large/15 by default?
lw=2;%1.5
tauend=size(A,2);
newend=min(tauend,500);
Ni=repmat(NN,1,tauend);
B=sum(Acum.*Ni,1)/sum(NN);
n=length(NN);

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
cmap=colormap(lines);
%cmap=.9*cmap;
colormap(cmap)
cc=colormap;
%imagesc(logN);
%cc=flipud(colormap);;
y1=Acum(sample,:); y2=A(sample,:); cc=[.447,.553,.647];
%
y1=y1(:,1:newend); y2=y2(:,1:newend);
T=1:newend;
hold on
%for i=1:thismany%Yes, this loop is clonky
    plot(T,y1,':','linewidth',lw)%,'color',cc(i,:));%[.5,.5,.5])%cc(i,:)); 'o-'
    plot(T,y2,'-','linewidth',lw)%,'color',cc(i,:));%[0,0,0])%cc(i,:)); 'o-'
%end
xlabel('Time (years)','FontSize',fs)
ylabel('Attack rate')%('Proportion immune','FontSize',fs) Relative attack
set(gca,'FontSize',fs);
maxY=max(max([y1;y2])); maxY=min(maxY+.1,1);%+.1
axis([1,newend,0,maxY])
grid on
grid minor
if n>1 && cbaroff~=1
    colorbar
    caxis([0,max(logN)])
end
hold off
%}
