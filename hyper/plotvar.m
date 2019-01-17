function plotvar(X,avg,R,reg,name)
figure;
hold on
N = length(X);
Xuniq = unique(X);
Xs = min(X):(max(X)-min(X))/9:max(X);

if length(Xuniq) <= 10
    boxplot(R,X,'positions',X+avg,'plotstyle','compact','jitter',0.3);
else
    scatter(X,R+normrnd(0,0.01,N,1))
end

%b = [ones(N,1) X]\R;
plot(Xs, [ones(size(Xs,2),1) (Xs-avg)']*reg);
title(name);
xlim([min(X) max(X)]);
hold off;
end

