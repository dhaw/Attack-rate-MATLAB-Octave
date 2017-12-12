function f=fghist(f,g)
l=length(f);
burn=10;
F=f(burn+1:end); G=g(burn+1:end);
maxz=max(f);
int=.01;
fs=15;
col1=[0,0,0]; col2=[.5,.5,.5];
figure
subplot(2,1,1)
[counts,centers]=hist(F,(0:int:maxz));
bar(centers,counts,'facecolor',col1,'edgecolor',col1,'barwidth',1);
%maxF=max(hist(F,(-int/2:int:maxz+int/2)));
maxF=max(counts);
axis([-int/2,maxz+int/2,0,maxF]);
%xlabel('Proportion immune','FontSize',fs)
ylabel('Frequency (total)','FontSize',fs)
subplot(2,1,2)
[counts1,centers1]=hist(G,(0:int:maxz));
bar(centers1,counts1,'facecolor',col2,'edgecolor',col2,'barwidth',1);
axis([-int/2,maxz+int/2,0,maxF]);
xlabel('Proportion immune','FontSize',fs)
ylabel('Frequency (new)','FontSize',fs)