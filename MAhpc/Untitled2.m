function f=plotMAdelay(F,G,fp)

[n,lt,thismany]=size(F);
Fcorr=nan(thismany,lt);
Gcorr=Fcorr;
Fvar=nan(thismany,1);
for i=1:lt
    parfor j=1:thismany
        ccf=corrcoef(F(:,fp