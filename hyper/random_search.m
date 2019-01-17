X_4C_8x10
opt = layoutopt(targets,weights);

%% run search
R = [];
D = [];

O = [];

for i=1:200
    o = randomopt(opt);
    ix = randperm(size(X0,1));
    X1 = X0(ix(1:o.population),:);
    [~,~,~,T,stats] = banditgasimulation(o,X1);
    %[~,~,~,T,stats] = batchgasimulation(o,X1);
    R = [R stats.expected_rewards(end)];
    D = [D stats.duplicates(end)];
    O = [O; o];
    i
end
%% save
results = struct();
results.R = R';
results.duplicates = D';
results.population = [O.population]';
results.mutation = [O.mutation]';
results.elitism = [O.elitism]';
results.alpha = [O.alpha]';
results.evaluation = {O.evaluation}';
results.crossover_rate = [O.crossover_rate]';
results.repetitions = [O.repetitions]';
results.exploration = [O.exploration]';
results.crossover_type = {O.crossover_type}';
results.creation = [O.creation]';

writetable(struct2table(results),'random_search.csv')