function f=plotMAdelay(F,G,fp)

[n,lt,thismany]=size(F);
Fcorr=nan(thismany,lt);
Gcorr=Fcorr;
tempv=zeros(thismany,2);
for i=1:lt
    tempvi=tempv;
    parfor j=1:thismany
        ccf=corrcoef(F(:,i,j),fp);
        tempv(j,1)=ccf(2);
        ccg=corrcoef(G(:,i,j),fp);
        tempvi(j,2)=ccg(2);
    end
    Fcorr(:,i)=tempvi(:,1);
    Gcorr(:,i)=tempvi(:,2);
end
Fvar=nanvar(Fcorr,[],2);
Gvar=nanvar(Gcorr,[],2);

fs=12; lw=1;%30;%Font size
figure
T=1:lt;
subplot(2,1,1)
hold on
plot(T,Fcorr,'-','linewidth',lw)%,'color',[.5,0,0])%'o-'
hold off
%xlabel('Time (years)','FontSize',fs)
ylabel('cc_{pand} (H1N1)','FontSize',fs)
%axis ([tmin-1,tmax,0,maxAttack]);
axis tight
set(gca,'FontSize',fs);
grid on
grid minor
box on
subplot(2,1,2)
hold on
plot(T,Gcorr,'-','linewidth',lw)%,'color',[0,0,.5])%'o-'
hold off
xlabel('Time (years)','FontSize',fs)
ylabel('cc_{pand} (H3N2)','FontSize',fs)
%axis ([tmin-1,tmax,0,maxAttack]);
axis tight
set(gca,'FontSize',fs);
grid on
grid minor
box on