function err = archRootFun(alpha,type,targetTau)
% if abs(alpha) < realmin
%     tau = 0;
% else
%     tau = 1 + 4 .* quadg(@lambdaarch,0,1,[],[],type,alpha);
% end

err = tau - targetTau;
