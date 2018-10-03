function f=sensAnalPlot(A,NN)
v=(.2:.2:4);
real=2.72;

A(A<0)=NaN;
rangeA=range(A,1);
NNrep=repmat(NN,1,size(A,2)); A(NNrep==0)=NaN;
repN=repmat(NN,1,size(A,2));
meanA=nansum((A.*repN),1)/sum(NN);
labels=cellstr(num2str(v'));
maxA=max(max(A));
maxA=max(1,maxA);

col1=[.165,.31,.431];%[1,0,0];%[0.0512,0.4600,0.8633];%[.165,.31,.431];
col2=[.447,.553,.647];
col3=[0,0,0];
col4=[.071,.208,.322];

figure
fs=18;
hold on
%plot(.52*ones(1,101),(0:maxA/100:maxA),'k.','linewidth',2)
plot([real,real],[0,maxA],'k--','linewidth',1.5)%
plot(v,meanA','-','color',col3,'linewidth',1.5);
bh=boxplot(A,'labels',labels,'position',v,'colors',col1,'symbol','.');%[.165,.31,.431][.447,.553,.647]
set(gca,'xtickmode','auto','xticklabelmode','auto')
axis([v(1)-.1,v(end)+.1,0,1]);%maxA
set(bh,'linewidth',1.5);
grid on
grid minor
xlabel('Distance power p','FontSize',fs); ylabel('Attack rate','FontSize',fs); set(gca,'FontSize',fs)%Proportion immune M_{80}
%Offset a Population power \alpha Distance power p