function f=plotBoxAgg(A,D)%(f,NN,D)%(B1,B2,B3,B4)%,NN
[a,b]=size(D);
v=(2:2:20); lv=length(v);%2:2:20
meanA=zeros(1,lv);
aggBefore=0;

if aggBefore==1
	D=B1+B2+B3+B4;
	c1=aggCells(B1,v(1)); [a,b]=size(c1);
	flist=NaN(a*b,lv);
else
    Abefore=A;
    f=A(:,1);
    flist=NaN(length(f),lv);
end

if aggBefore==1
    for i=1:lv
        c1=aggCells(B1,v(i));
        c2=aggCells(B2,v(i));
        c3=aggCells(B3,v(i));
        c4=aggCells(B4,v(i));
        [gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,ages0]=prepFluAge(c1,c2,c3,c4,1.8,0,v(i)/10);
        [f,g]=finalSizeMulti(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaD,1,2,10^(-16));
        flist(1:length(f),i)=f;%.*NN;%********test case
        meanA(i)=nanmean(f.*NN);
    end
else
    %c1=aggCells(B1,v(1));
    %c2=aggCells(B2,v(1));
    %c3=aggCells(B3,v(1));
    %c4=aggCells(B4,v(1));
    %[gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,ages0]=prepFluAge(c1,c2,c3,c4,1.8,0,v(1));
    %[f,g]=finalSizeMulti(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaD,1,2,10^(-16));
    flist(1:length(f),1)=f;
    for i=2:lv
        [fnew,NNnew]=aggResult(D,f,v(i)/2);
        flist(1:length(fnew),i)=fnew;
        meanA(i)=nanmean(fnew.*NNnew);
    end
end
%A=padcat(v);
A=flist;
f=A;
%}

%rangeA=range(A,1);
%repN=repmat(NN,1,size(A,2));
%meanA=nanmean(A.*repN,1);
meanA=nanmean(A,1);

quantA=quantile(A,[0,.25,.5,.75,1],1);

labels=cellstr(num2str(v'));
maxA=max(max(A))+.1;
maxA=min(1,maxA);

fs=12; lw=1;
%{
figure
hold on
plot(v,meanA','-','color',[0,.25,.25],'linewidth',lw)
%plot(v,quantA,'o-','linewidth',lw)
fanChart(v,quantA');
%boxplot(A,'labels',labels,'position',v,'colors',[0.0512,0.4600,0.8633],'PlotStyle','compact','Notch','off','whisker',100000);%plot(v,meanA,'-','color',[0,.5,.5],'linewidth',2) %'labels',num2str(v),'position',v,

axis([v(1)-2,v(end)+2,0,min(maxA+.1,1)]);%maxA
grid on
grid minor
xlabel('Resolution (x100m)','FontSize',10); ylabel('Attack rate','FontSize',fs); set(gca,'FontSize',10)%Proportion immune M_{80}
%}
figure
maxA=max(max([Abefore;A]));
maxA=min(maxA+.1,1);
minA=min(min([Abefore;A]));
minA=max(minA-.1,0);
subplot(2,1,1)
hold on
fanChart(v,Abefore');
%plot(v,meanA','-','color',[0,.25,.25],'linewidth',lw)
hold off
axis([v(1)-2,v(end)+2,minA,maxA]);
grid on
grid minor
title('Aggregation before simulation')  
xlabel('Resolution (x100m)'); ylabel('Attack rate');
set(gca,'FontSize',fs)
subplot(2,1,2)
hold on
fanChart(v,A');
%plot(v,meanA','-','color',[0,.25,.25],'linewidth',lw)
hold off
axis([v(1)-2,v(end)+2,minA,maxA]);
grid on
grid minor
title('Aggregation after simulation') 
xlabel('Resolution (x100m)'); ylabel('Attack rate');
set(gca,'FontSize',fs)