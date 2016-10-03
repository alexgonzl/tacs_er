function x2 = vecRej(x,y)
    xm = x-nanmean(x);
    ym = y-nanmean(y);
    x1=nansum(xm.*ym)/norm2(ym)^2*ym;
    x2=x-x1;
end
function x2=norm2(x)
    x2=sqrt(nansum(x.^2));
end