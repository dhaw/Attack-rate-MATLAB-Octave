function f=plotBoxDelta(B1,B2,B3,B4)%,NN %(A,D) (f,NN,D)
%[a,b]=size(D);
delta=(0:.05:1); ld=length(delta);%2:2:20
meanA=zeros(1,ld);
v=(2:2:6); lv=length(v);
meanA2=zeros(1,lv);

	c1=aggCells(B1,2); [a,b]=size(c1);
    c2=aggCells(B2,2);
    c3=aggCells(B3,2);
    c4=aggCells(B4,2);
	flist=NaN(a*b,ld);
    flist2=NaN(a*b,lv);
    

for i=1:ld
        [gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,ages0]=prepFluAge(c1,c2,c3,c4,1.8,0,delta(i));
        [f,g]=finalSizeMulti(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaD,1,2,10^(-16));
        flist(1:length(f),i)=f;%.*NN;%********test case
        meanA(i)=nanmean(f.*NN);
        if i==ld
            flist2(1:length(f),1)=f;
            meanA2(1)=nanmean(f.*NN);
        end
end
A=flist;
f=A;
%{
for i=2:lv
        c1=aggCells(B1,v(i));
        c2=aggCells(B2,v(i));
        c3=aggCells(B3,v(i));
        c4=aggCells(B4,v(i));
        [gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,ages0]=prepFluAge(c1,c2,c3,c4,1.8,0,v(i)/10);
        [f,g]=finalSizeMulti(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaD,1,2,10^(-16));
        flist2(1:length(f),i)=f;%.*NN;%********test case
        meanA2(i)=nanmean(f.*NN);
end
A2=flist2;
g=A2;
%}
%{
meanA=nanmean(A,1);
quantA=quantile(A,[0,.25,.5,.75,1],1);
labels=cellstr(num2str(delta'));
maxA=max(max(A))+.1;
maxA=min(1,maxA);
%}
fs=12; lw=1;
figure
%subplot(2,1,1)
maxA=max(max(A));
maxA=min(maxA+.1,1);
minA=min(min(A));
minA=max(minA-.1,0);
hold on
fanChart(delta,A');
%plot(v,meanA','-','color',[0,.25,.25],'linewidth',lw)
hold off
axis([0,1,minA,maxA]);
grid on
grid minor
title('Limited mobility')  
xlabel('\delta'); ylabel('Attack rate');
set(gca,'FontSize',fs)
%{
subplot(2,1,2)
hold on
fanChart(v,A2');
plot(v,meanA2','-','color',[0,0,0],'linewidth',lw)
hold off
axis([0,v(end),minA,maxA]);
grid on
grid minor
title('Aggregation before simulation')  
xlabel('Resolution (x100m)'); ylabel('Attack rate');
set(gca,'FontSize',fs)
%}