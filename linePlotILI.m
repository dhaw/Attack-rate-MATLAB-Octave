<<<<<<< HEAD
function f=linePlotILI(data)
A=table2array(data(:,2:4));
remove=sum(isnan(A),2);
index=find(remove==0);
A=A(index,2:end);

fs=15;%30;%Font size
tauend=size(A,2);
%
figure
T=1:tauend;
hold on
colormap grey
plot(T,A,'o:','linewidth',1.5,'color','k')
xlabel('Time (years)','FontSize',fs)
ylabel('Proportion immune','FontSize',fs)
set(gca,'FontSize',fs);
maxY=max(max([Acum;A])); maxY=min(maxY+.1,1);
axis([1,tauend,0,maxY])
grid on
=======
function f=linePlotILI(data)
A=table2array(data(:,2:4));
remove=sum(isnan(A),2);
index=find(remove==0);
A=A(index,2:end);

fs=15;%30;%Font size
tauend=size(A,2);
%
figure
T=1:tauend;
hold on
colormap grey
plot(T,A,'o:','linewidth',1.5,'color','k')
xlabel('Time (years)','FontSize',fs)
ylabel('Proportion immune','FontSize',fs)
set(gca,'FontSize',fs);
maxY=max(max([Acum;A])); maxY=min(maxY+.1,1);
axis([1,tauend,0,maxY])
grid on
>>>>>>> acfad3e4da850e525286e08138beef9f06a8c2cc
hold off