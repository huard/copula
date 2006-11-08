function err = frankRootFun(alpha,targetTau)
if abs(alpha) < realmin
    tau = 0;
else
    tau = 1 + 4 .* (debye1(alpha)-1) ./ alpha;
end
err = tau - targetTau;