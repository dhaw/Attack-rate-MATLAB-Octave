function f=RPvarPrev(tvec,Y,boxLat,boxLong)
%Y(Y==0)=NaN;
Yvar=nanvar(Y,0,2);
Yvar=sqrt(Yvar);
fs=12; lw=2;
logVar=1;
logDist=1;
distNotSpeed=1;%1 for distance to peak, 0 (or other) for speed
incidence=0;
%
plotMax=max(Yvar);
figure
if logVar==1
    semilogy(tvec,Yvar,'r-','linewidth',lw)
    plotMin=min(Yvar(Yvar~=0));
    logPlotMax=ceil(log10(plotMax));
    logPlotMin=floor(log10(plotMin));
    tickvec=logspace(logPlotMin,logPlotMax,logPlotMax-logPlotMin+1);
    plotMax=10^(logPlotMax);
    axis([0,tvec(end),plotMin,plotMax])
    xlabel('Time (days)','FontSize',fs);
    if incidence==1
        ylabel('\sigma(inc.)','FontSize',fs);
    else
        ylabel('\sigma(prev.)','FontSize',fs);
    end
    set(gca,'FontSize',fs,'ytick',tickvec);
    grid on
    grid minor
else
    plot(tvec,Yvar,'r-','linewidth',lw)
    axis ([0,tvec(end),0,plotMax])
    xlabel('Time (days)','FontSize',fs);
    if incidence==1
        ylabel('\sigma(inc.)','FontSize',fs);
    else
        ylabel('\sigma(prev.)','FontSize',fs);
    end
    set(gca,'FontSize',fs);
    grid on
    grid minor
end
%}
%
n=boxLat*boxLong;
nums=reshape((1:n),boxLat,boxLong);
line=flipud(nums(:,1));
%Line to plot - box-specific
Ysub=Y(:,line);
lt=length(tvec);
peakLoc=zeros(lt,1);
for i=1:lt
    yi=Ysub(i,:);
    [maxi,indi]=max(yi);
    peakLoc(i)=indi;
end

if distNotSpeed~=1
    peakLoc=diff(peakLoc).*diff(tvec); tvec(1)=[];
end

plotMax=max(peakLoc);
figure
if logDist==1
    semilogy(tvec,peakLoc,'r-','linewidth',lw)
    logPlotMin=-1; plotMin=10^logPlotMin;
    logPlotMax=ceil(log10(plotMax));
    tickvec=logspace(logPlotMin,logPlotMax,logPlotMax-logPlotMin+1);
    plotMax=10^(logPlotMax);
    axis([0,tvec(end),plotMin,plotMax])
    xlabel('Time (days)','FontSize',fs);
    if distNotSpeed==1
        ylabel('Distance to peak','FontSize',fs);
    else
        ylabel('Peak speed','FontSize',fs);
    end
    set(gca,'FontSize',fs,'ytick',tickvec);
    grid on
    grid minor
else
    plot(tvec,peakLoc,'r-','linewidth',lw)
    axis([0,tvec(end),0,plotMax])
    xlabel('Time (days)','FontSize',fs);
    if distNotSpeed==1
        ylabel('Distance to peak','FontSize',fs);
    else
        ylabel('Peak speed','FontSize',fs);
    end
    set(gca,'FontSize',fs);
    grid on
    grid minor
end
%}