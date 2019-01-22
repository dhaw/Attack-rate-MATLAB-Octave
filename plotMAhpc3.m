function f=plotMAhpc3(X)

[a,b,c]=size(X);

%Concatenate simulations:
%X=reshape(X,a,b*c,1);

%Distribution of variances:
X=nanvar(X,[],2);
X=reshape(X,[a,c,1]);

eps=(0:.1:1);
leps=length(eps);
xbox=prctile(X,[0,25,50,75,100],2);
%ybox=prctile(Y,[0,25,50,75,100],2);
%}
fs=12; lw=2;
figure;
%hold on
%plot(1:leps,y,'k-','linewidth',lw)
%plot(1:leps,z,'k--','linewidth',lw)%,'color',[0,.5,0]
h=boxplot(xbox','whisker',10);
set(h,{'linew'},{1})
set(gca,'xtick',(1:2:leps),'xticklabel',eps(1:2:leps));
hold off
set(gca,'fontsize',fs)
axis ([0,leps+1,-1,1])
xlabel('\epsilon')
ylabel('Correlation')
grid on
grid minor
box on
%}
%%Drift - rank order and epsilon:
%

%}