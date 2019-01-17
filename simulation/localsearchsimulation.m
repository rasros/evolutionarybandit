function [X,Z,Ze,T,stats,exitflag] = localsearchsimulation(opt)

modtick = max(1,ceil(opt.T/100));
stats = struct();
stats.mean_score = zeros(1,modtick);
stats.rewards = zeros(1,modtick);
stats.max_score = zeros(1,modtick);

x0 = opt.generate();
evals = zeros(1,modtick);
tevals = 0;
T = (1:modtick)*opt.population * opt.evaluations;

    function f = fit(x)
        tevals = tevals + opt.evaluations;
        t = ceil(tevals / modtick);
        evals(t) = evals(t) + opt.evaluations;
        s = opt.fitness(round(x));
        r = opt.experiment(s,opt.evaluations);
        f = 1 - r/opt.evaluations;
        stats.mean_score(t) = stats.mean_score(t)+s*opt.evaluations;
        stats.rewards(t) = stats.rewards(t)+r;
        stats.max_score(t) = max(stats.max_score(t), s);
    end

paoptions = optimoptions('patternsearch', 'TolMesh', 0.9, 'ScaleMesh', false, ...
    'MaxTime', inf, 'MaxFunctionEvaluations', inf, 'SearchFcn', []);

[X,Ze,exitflag] = patternsearch(@fit,x0,[],[],[],[],1,opt.C,[],paoptions);
X = round(X);
Z = opt.fitness(X);

t1 = ceil(tevals / modtick);

stats.duplicates = zeros(1,t1);
stats.mutations = ones(1,t1);
stats.mutp = ones(1,t1);
stats.pop_size = ones(1,t1);

stats.mean_score = stats.mean_score ./ evals;
stats.mean_reward = stats.rewards ./ evals;

stats.mean_score(t1+1:end) = [];
stats.mean_reward(t1+1:end) = [];
stats.rewards(t1+1:end) = [];
stats.max_score(t1+1:end) = [];
T(t1+1:end) = [];

a = opt.prior_succ;
b = opt.prior_trial;
Ze = Ze*opt.evaluations+a;
Ze = Ze/(b+opt.evaluations);

end