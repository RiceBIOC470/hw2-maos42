function mm=meannonan(x)
%notin=isnan(x) | isinf(x);
%x(notin)=[];
n=length(x);
mm=nanmean(x(1:n,1:n));

