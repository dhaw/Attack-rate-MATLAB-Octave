%%
%Agg before
%
function f=aggBox(C1,C2,C3,C4,fc)%(C1,C2,C3,C4,fc)
fcCelldist=.2;
fc1=fc;
aggfact=(1:20);
lagg=length(aggfact);
X=zeros(5,lagg);
X(:,1)=prctile(fc,[0,25,50,75,100]);
for i=2:lagg
    aggi=aggfact(i);
    celldist=fcCelldist*aggi;
    c1=aggCells(C1,aggi);
    c2=aggCells(C2,aggi);
    c3=aggCells(C3,aggi);
    c4=aggCells(C4,aggi);
    %Call=C1+C2+C3+C4;
    beta3=1;
    [gamma,NN,n,nbar,na,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,ages0]=prepFluAge(c1,c2,c3,c4,1.8,0,1,celldist,13.5);
    [fc,~,~]=finalSizeMulti(gamma,n,nbar,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3,1,2,.00001,0);
    xi=prctile(fc,[0,25,50,75,100]);
    X(:,i)=xi;
end
f=X;
figure;
%plot(0:6,HH(:,1),'k-','linewidth',1)
hold on
h=boxplot(X,'whisker',10);%,'positions',res,'labels',res);
set(h,{'linew'},{1})
set(gca,'xtick',(2:2:lagg),'xticklabel',fcCelldist*(2:2:lagg));
hold off
set(gca,'fontsize',12)
%set(gca,'xtick',0:6)
axis ([aggfact(1)-.5,aggfact(end)+.5,min(min(X))-.05,max(max(X))+.05])%max(min(fc1)-.05,0),max(fc1)+.05])
title('Before')
xlabel('Aggregation (km)')
ylabel('Attack rate')
grid on
grid minor
box on
%}
%%
%Agg after
%{
function f=aggBox(fc,Call)
fcCelldist=.2;
res=(1:20);
lres=length(res);
X=zeros(5,lres);
X(:,1)=prctile(fc,[0,25,50,75,100]);
for i=2:lres
    [a,~]=aggResult(Call,fc,res(i));
    xi=prctile(a,[0,25,50,75,100]);
    X(:,i)=xi;
end

figure;
%plot(0:6,HH(:,1),'k-','linewidth',1)
hold on
h=boxplot(X),'whisker',10);%,'positions',res,'labels',res);
set(h,{'linew'},{1})
set(gca,'xtick',(2:2:lres),'xticklabel',fcCelldist*(2:2:lres));
hold off
set(gca,'fontsize',12)
%set(gca,'xtick',0:6)
axis ([res(1)-.5-fcCelldist,res(end)+.5,min(min(X))-.05,max(max(X))+.05])%max(min(fc)-.05,0),max(fc)+.05])
title('After')
xlabel('Aggregation (km)')
ylabel('Attack rate')
grid on
grid minor
box on
%}