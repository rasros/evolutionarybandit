function barvar(X,R,name)
figure
boxplot(R,X,'plotstyle','compact','jitter',0.3);
title(name);
end