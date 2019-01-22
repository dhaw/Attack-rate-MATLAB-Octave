function f=plotMAhpc2(X,Y,Z)
%%Theta geom:
%
eps=(0:.1:1);
leps=length(eps);
x=X(:,:,1);
xbox=prctile(x,[0,25,50,75,100],1);

y=min(Y,[],2);
%z=Z(:,1);%min(Z,[],2);
z=var(X,[],3);
z=max(z,[],1);

fs=12; lw=2;
figure;
hold on
plot(1:leps,y,'k-','linewidth',lw)
plot(1:leps,z,'k--','linewidth',lw)%,'color',[0,.5,0]
h=boxplot(xbox,'whisker',10);
set(h,{'linew'},{1})
set(gca,'xtick',(1:2:leps),'xticklabel',eps(1:2:leps));
hold off
set(gca,'fontsize',fs)
axis ([0,leps+1,0,1])
ylabel('Attack rate')
grid on
grid minor
box on
%}
%%Drift - rank order and epsilon:
%
%}