function [OUT] = TDNL_aux(U, KK, c, PP)

%  negels
nU = U(U<0);
nei = find(U<0);

if (min(nei)==1)
    innern = nU.^2 - U(nei).^2;
    Xn = nU + c.*(innern)./PP;
    Xn(1) = U(1) + c*(U(1)*U(1) - U(end)*U(end))/PP;
else
    innern = nU.^2 - U(nei-1).^2;
    Xn = nU + c.*(innern)./PP;
end


% posels
pU = U(U>0);
poi = find(U>0);

if (max(poi)==length(U))
    innerp = (U(poi).^2) - pU.^2;
    Xp = pU + c*(innerp)/PP;
    Xp(end) = U(end) + c*(U(1)*U(1) - U(end)*U(end))/PP;
else
    innerp = (U(poi+1).^2) - pU.^2;
    Xp = pU + c*(innerp)/PP;
end


OUT = zeros(size(U));
OUT(nei) = Xn;
OUT(poi) = Xp;

end