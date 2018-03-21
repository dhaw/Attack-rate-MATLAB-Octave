function f=sensAnalPlot(A,NN)
v=(.2:.2:4);
real=2.72;

A(A<0)=NaN;
rangeA=range(A,1);
NNrep=repmat(NN,1,size(A,2)); A(NNrep==0)=NaN;
repN=repmat(NN,1,size(A,2));
meanA=nansum(A.*repN,1)/sum(NN);
labels=cellstr(num2str(v'));
maxA=max(max(A));
maxA=max(1,maxA);
figure
fs=20;
hold on
%plot(.52*ones(1,101),(0:maxA/100:maxA),'k.','linewidth',2)
plot([real,real],[0,maxA],'k--','linewidth',1.5)%
plot(v,meanA','-','color',[.071,.208,.322],'linewidth',1.5);
bh=boxplot(A,'labels',labels,'position',v,'colors',[.447,.553,.647],'symbol','.');%[.165,.31,.431][.447,.553,.647]
set(gca,'xtickmode','auto','xticklabelmode','auto')
axis([0,v(end),0,1]);%maxA
set(bh,'linewidth',1.5);
grid on
xlabel('Distance power p','FontSize',fs); ylabel('Proportion immune','FontSize',fs); set(gca,'FontSize',fs)%Proportion immune M_{80}
%Offset a Population power \alpha Distance power p