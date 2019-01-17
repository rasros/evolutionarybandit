function [X,Z,Ze,T,stats,exitflag,output] = batchgasimulation(opt,X0)
generations = ceil(opt.T / (opt.population * opt.repetitions));
stats = struct();
stats.mean_score = zeros(1,generations);
stats.rewards = zeros(1,generations);
stats.max_score = zeros(1,generations);
stats.duplicates = zeros(1,generations);
stats.evictions = ones(1,generations);
T = (1:generations)*opt.population*opt.repetitions;
T(end) = opt.T;
totalF = 0;
if nargin < 2
    X = zeros(opt.population,opt.N);
    for i=1:opt.population
        X(i,:) = opt.generate();
    end
else
    X = X0(1:opt.population,:);
end

gen = 1;

    function f = fit(x)
        totalF = min(totalF + opt.repetitions,opt.T);
        if totalF < opt.T
            s = opt.fitness(x);
            r = opt.experiment(s,opt.repetitions);
            f = 1 - r/opt.repetitions;
            stats.mean_score(gen) = stats.mean_score(gen)+s;
            stats.rewards(gen) = stats.rewards(gen)+r;
            stats.max_score(gen) = max(stats.max_score(gen), s);
        else
            f = nan;
        end
    end

    function children = mutation(opt, parents, s, X)
        stats.duplicates(gen) = 1-(size(unique(X,'rows'),1) / size(X,1));
        children = zeros(size(parents,1),size(parents,2));
        for i=1:size(parents,1)
            if rand() < opt.mutation / (size(parents,1)/opt.population)
                children(i,:) = opt.mutate(parents(i,:));
            else
                children(i,:) = parents(i,:);
            end
        end
        gen = s.Generation+1;
    end

gaoptions = optimoptions('ga',...
    'MaxGenerations', max(1,generations-1), ...
    'MaxStallGenerations', Inf, ...
    'MaxTime', Inf, ...
    'FunctionTolerance', 0, ...
    'PopulationSize', opt.population, ...
    'CreationFcn', @(a1,a2,a3) X, ...
    'EliteCount', min(opt.population-1,opt.elitism), ...
    'CrossoverFraction', opt.crossover_rate, ...
    'CrossoverFcn', @(parents,a1,a2,a3,a4,X) crossover(parents,opt,X), ...
    'MutationFcn', @(parents,a1,a2,a3,o,scores,a6) mutation(opt,a6(parents,:),o,X));

[~,~,exitflag,output,X,Ze] = ga(@fit,opt.N,[],[],[],[],1,opt.C,[],[],gaoptions);

Z = sortrows([Ze (1:length(Ze))'], 1);
X = X(Z(:,2),:);
Ze = (1-Ze(Z(:,2)))*opt.repetitions;

Z = zeros(size(X,1),1);
for i=1:size(X,1)
    Z(i) = opt.fitness(X(i,:));
end
stats.mean_score = stats.mean_score / opt.population;
stats.mean_score(end) = mean(Z);
stats.duplicates(end) = 1-size(unique(X,'rows'),1) / size(X,1);
stats.mean_reward = stats.rewards ./ [T(1) diff(T)];
stats.expected_rewards = stats.mean_score;
end

function Xp = crossover(parents, lopt, X)
nk = length(parents)/2;
Xp = zeros(nk,lopt.N);
index = 1;
for i=1:nk
    x1 = X(parents(index),:);
    x2 = X(parents(index+1),:);
    Xp(i,:) = lopt.crossover(x1,x2);
    index = index+2;
end
end