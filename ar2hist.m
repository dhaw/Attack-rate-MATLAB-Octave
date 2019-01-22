function f=ar2hist(fc,NN)
meanfc=sum(fc.*NN)/sum(NN);
fs=12;
lw=2;
figure
hist(fc,20,'facecolor',[0,0,.5]);
hold on
plot([meanfc,meanfc],[0,1000],'--','color',[.5,0,0],'linewidth',lw)

hold off
xlabel('Attack rate')
ylabel('Frequency')
set(gca,'fontsize',fs)
axis([min(fc)-.01,max(fc)+.01,0,max(hist(fc,20))])
grid on
grid minor
box on
%f=h1;