function opt2 = randomopt(opt)
opt2=opt;
maxpop = min(ceil(opt.T/5), 99);
opt2.population = 1+randi(maxpop);
opt2.mutation = rand();
opt2.repetitions = randi(floor(opt2.T/opt2.population));
opt2.elitism = randi(opt2.population)-1;
opt2.crossover_rate = rand();
opt2.alpha = rand();
opt2.creation = rand()*0.05;

evaluations = {'thompson','thompson_optimistic','ucb'};
opt2.evaluation = evaluations{randi(length(evaluations))};
opt2.exploration = rand()*20;

crosstypes = {'scatter','single_point','rectangle'};
opt2.crossover_type = crosstypes{randi(length(crosstypes))};
end
