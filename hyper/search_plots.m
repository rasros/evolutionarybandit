batch = 0;

%% load
if batch
    T = readtable('random_search_batch.csv');
else
    T = readtable('random_search_bandit.csv');
end

N = length(T.R);

TS = sortrows([(1:size(T,1))' T.R], -2);
TS = TS(1:size(T,1),1);
T = T(TS,:);

%% regression model
% Models are found through factor analysis and deduction
if batch
    
    % Main effects
    %   population, mutation, elitism, crossover_rate, repetitions,
    % Interaction effects
    %   population * repetitions
    %   repetitions * crossover_rate
    %   population * mutation * elitism
    
    % population = 20
    % mutation = 0.1
    % elitism = 5
    % crossover_rate = 0.8
    % repetitions = 10
        
    X1 = [log10(T.population) T.mutation log10(1+T.elitism) T.crossover_rate log10(T.repetitions)];
    Xmean = mean(X1);
    X1 = bsxfun(@minus,X1,Xmean);
    X2 = [X1(:,1).*X1(:,5) ...
        X1(:,4).*X1(:,5) ...
        X1(:,1).*X1(:,2).*X1(:,3)];
    X = [ones(size(T,1),1) X1 X2];
else
    % Main effects
    %   population, mutation, elitism, alpha, exploration, creation
    % Interaction effects
    %   population * elitism
    %   population * alpha
    %ix = arrayfun(@strcmp,T.evaluation,repmat({'thompson'},size(T,1),1));
    %T = T(ix,:);
    
    X1 = [log10(T.population) T.mutation log10(1+T.elitism) T.alpha T.exploration T.creation];
    Xmean = mean(X1);
    X1 = bsxfun(@minus,X1,Xmean);
    
    X2 = [X1(:,1).*X1(:,3) X1(:,1).*X1(:,4)];
    X = [ones(size(T,1),1) X1 X2];
end

model = regress(T.R,X);

%% plots

plotvar(log10(T.population),Xmean(1),T.R,[model(1); model(2)],'population');
plotvar(T.mutation,Xmean(2),T.R,[model(1); model(3)],'mutation');
plotvar(log10(1+T.elitism),Xmean(3),T.R,[model(1); model(4)],'elitism');
barvar(T.crossover_type,T.R,'crossover type');

if batch
    plotvar(log10(T.repetitions),Xmean(5),T.R,[model(1); model(6)],'repetitions');
    plotvar(T.crossover_rate,Xmean(4),T.R,[model(1); model(5)],'crossover rate');
    plotvar(X2(:,1),0,T.R,[model(1); model(7)],'population*repetitions');
    plotvar(X2(:,2),0,T.R,[model(1); model(8)],'repetitions*crossover rate');
    plotvar(X2(:,3),0,T.R,[model(1); model(9)],'population*mutation*elitism');
else
    plotvar(T.alpha,Xmean(4),T.R,[model(1); model(5)],'alpha');
    plotvar(T.exploration,Xmean(5),T.R,[model(1); model(6)],'exploration');
    plotvar(T.creation,Xmean(6),T.R,[model(1); model(7)],'creation');
    
    plotvar(X2(:,1),0,T.R,[model(1); model(8)],'population*elitism');
    plotvar(X2(:,2),0,T.R,[model(1); model(9)],'population*alpha');
    barvar(T.evaluation,T.R,'evaluation');
end