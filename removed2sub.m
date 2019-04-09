%
        lmax=max(maxYears1,maxYears2);
        X1=X; Y1=Y;
        %X=X.*Nratio1; Y=Y.*Nratio2;
        if maxYears1<lmax
            X1(:,end+1:lmax)=0;
            Nratio1=Nratio2;
        elseif maxYears2<lmax
            Y1(:,end+1:lmax)=0;
            Nratio2=Nratio1;
        end
        Sdiff=X1+Y1-1; Sdiff(Sdiff<0)=0;
        findSdiff=find(Sdiff); ldiff=length(findSdiff);
        randNZ2=rand(length(findSdiff),1);
        
        %sub1=Sdiff; sub2=Sdiff;
        %sub1(findSdiff)=sub1(findSdiff).*randNZ2;
        %sub2(findSdiff)=sub2(findSdiff).*(1-randNZ2);
        
        Xtake=X1(findSdiff); %sub1(Sdiff==0)=0; %Xtake is TOTAL X in relevant cells
        Ytake=Y1(findSdiff); %sub2(Sdiff==0)=0;
        denom=Xtake+Ytake;
        w1=Xtake./denom.*Sdiff(findSdiff); w2=Ytake./denom.*Sdiff(findSdiff);
        %r1=rand(ldiff,1); r2=rand(ldiff,1);
        %sub1=Xtake.*(w1.*(1-r1)+r1); sub2=Ytake.*(w2.*(1-r2)+r2); 
        sub1=Xtake.*w1; sub2=Ytake.*w2; 
        X1(findSdiff)=X1(findSdiff)-sub1; Y1(findSdiff)=Y1(findSdiff)-sub2;
        %
        %X1=X1-sub1; Y1=Y1-sub2; 
        
        X=X1(:,1:maxYears1); Y=Y1(:,1:maxYears2);
        %b1=[znbar,X,zeros(nbar,lp1-maxYears1)]; b2=[Y,zeros(nbar,lp2-maxYears2)];
        b1=[znbar,X]; b2=[znbar,Y];
        %b1(isinf(b1)==1)=0; b1(isnan(b1)==1)=0;
        %b2(isinf(b2)==1)=0; b2(isnan(b2)==1)=0;
        %
        %b1=[zeros(nbar,1),b1];
        %b2=[zeros(nbar,1),b2];
        X1=X1.*Nratio1; %X1(isinf(X1)==1)=0; X1(isnan(X1)==1)=0;
        Y1=Y1.*Nratio2; %Y1(isinf(Y1)==1)=0; Y1(isnan(Y1)==1)=0;
        Z1=zeros(nbar,lmax); Z1(findSdiff)=sub1+sub2;
        %Z1=sub1+sub2;
        if maxYears1<lmax
            Nratio3=Nratio1;
        else
            Nratio3=Nratio2;
        end
        Z1=Z1.*Nratio3;
        s01=sum(X1,2).*NNbar; s10=sum(Y1,2).*NNbar; s00=sum(Z1,2).*NNbar; %s11=NNbar./NNrep-s01-s10-s00; s11(NNrep==0)=0;
        s11=NNbar-s01-s10-s00; s11(NNrep==0)=0;