function err = archRootFun(alpha,type,targetTau)
if abs(alpha) < realmin
    tau = 0;
else
    tau = 1 + 4 .* quadl(@lambdaarch,0,1,[],[],alpha,type);
end
err = tau - targetTau;