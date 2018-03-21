function f=fghist(f,g)
l=length(f);
burn=10;
F=f(burn+1:end); G=g(burn+1:end);
maxz=max(f);
int=.01;
fs=15;
figure
subplot(2,1,1)
hist(F,(0:int:maxz));
maxF=max(hist(F,(-int/2:int:maxz+int/2)));
axis([-int/2,maxz+int/2,0,maxF]);
%xlabel('Proportion immune','FontSize',fs)
ylabel('Frequency (total)','FontSize',fs)
subplot(2,1,2)
hist(G,(0:int:maxz));
axis([-int/2,maxz+int/2,0,maxF]);
xlabel('Proportion immune','FontSize',fs)
ylabel('Frequency (new)','FontSize',fs)