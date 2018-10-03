function f=MAdata(X)
%Written quickly for ECMTB - some inpits need automating
years=(2009:2017);

fs=12; ms=5;%12; 20
lwx=1; lw=1.5; %lwx=1.5
%
%Several lines:
cmap=lines(6);
figure
plot(years,X,'o-','linewidth',lw,'markersize',ms)
xlabel('Year')
ylabel('Proportion')
set(gca,'fontsize',fs)
grid on
grid minor
legend('GBR','FRA','DEU','ESP','PRT','ITA','location','northeastoutside');
%}
%{
%Bars:
figure
Xsum=sum(X,2); X=X./repmat(Xsum,1,6);
bar(years,X,'stacked')
xlabel('Year')
ylabel('Proportion')
set(gca,'fontsize',fs)
grid on
grid minor
legend('H1','H3','Hna','BVic','BYam','Bna','location','northeastoutside')
%}