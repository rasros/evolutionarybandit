function [X,Z,Ze,T,stats] = banditgasimulation(opt,X0)

P = opt.population;
if nargin < 2
    X = zeros(P,opt.N);
    for i=1:P
        X(i,:) = opt.generate();
    end
else
    X = X0(1:P,:);
end
Z = cell(P,1);
statn = 20;
modtick = max(1,ceil(opt.T/statn));
T = unique([modtick:modtick:opt.T opt.T]);
statn = length(T);

stats = struct();
stats.expected_rewards = zeros(1,statn);
stats.mean_score = zeros(1,statn);
stats.rewards = zeros(1,statn);
stats.max_score = zeros(1,statn);
stats.duplicates = zeros(1,statn);
stats.mutations = zeros(1,statn);
stats.evictions = zeros(1,statn);
stats.threshold = zeros(1,statn);
stats.duplicate = zeros(1,statn);
TN = zeros(1,statn);

for i=1:P
    Z{i} = stat_struct(opt.fitness(X(i,:)));
    Z{i}.r = 1;
    Z{i}.n = 1;
end

for t=1:opt.T
    stati =  ceil(t / modtick);
    Zm = cell2mat(Z);
    TN(stati) = 1+TN(stati);
    
    %% Pick experiment genome
    bestV = -1;
    ix = -1;
    for i=1:size(X,1)
        z = Z{i};
        stats.max_score(stati) = max(z.f, stats.max_score(stati));
        v = evaluation(opt,opt.evaluation,z,Zm);
        if v > bestV
            bestV = v;
            ix = i;
        end
    end
    
    %% Perform experiment
    fn = opt.experiment(Z{ix}.f,1);
    stats.rewards(stati) = stats.rewards(stati)+fn;
    stats.expected_rewards(stati) = stats.expected_rewards(stati)+Z{ix}.f;
    if fn == 0 && Z{ix}.r == 0
        stats.evictions(stati) = stats.evictions(stati)+1;
        Z(ix) = [];
        Zm(ix,:) = [];
        X(ix,:) = [];
    else
        Z{ix}.r = Z{ix}.r+fn;
        Z{ix}.n = Z{ix}.n+1;
        [phat,ci]=binofit(Z{ix}.r,Z{ix}.n,opt.alpha);
        Z{ix}.phat = phat;
        Z{ix}.ci = ci;
    end
    
    %% Perform elimination
    CI = reshape([Zm.ci],2,size(X,1));
    threshold = max(CI(1,:));
    if opt.elitism > 0
        E = sortrows([Zm.r]',-1);
        eliteThreshold = E(min(size(X,1),opt.elitism));
    else
        eliteThreshold = inf;
    end
    
    stats.threshold(stati) = threshold;
    
    elim = find(CI(2,:) < threshold & [Zm.r] < eliteThreshold);
    
    if ~isempty(elim)
        stats.evictions(stati) = stats.evictions(stati)+length(elim);
        X(elim,:) = [];
        Z(elim) = [];
        Zm(elim) = [];
        CI(:,elim) = [];
    end
    
    %% Create new individuals
    P = max(opt.population,size(X,1));
    if rand() < opt.creation
        P = P+1;
    end
    X2 = zeros(P-size(X,1),opt.N);
    Z2 = cell(P-size(X,1),1);
    for k=1:P-size(X,1)
        if size(X,1) > 1
            ix1 = randi(size(X,1));
            tries1=1;
            while tries1<100 && CI(1,ix1) / threshold < rand()
                ix1 = randi(size(X,1));
                tries1=tries1+1;
            end
            ix2 = randi(size(X,1));
            tries1=1;
            while tries1<100 && ix1==ix2 || CI(1,ix2) / threshold < rand()
                ix2 = randi(size(X,1));
                tries1=tries1+1;
            end
            X2(k,:) = opt.crossover(X(ix1,:),X(ix2,:));
            if rand() < opt.mutation
                X2(k,:) = opt.mutate(X2(k,:));
                stats.mutations(stati) = stats.mutations(stati)+1;
            end
        else
            X2(k,:) = opt.mutate(X(1,:));
        end
        Z2{k} = stat_struct(opt.fitness(X2(k,:)));
    end
    X = [X; X2];
    Z = [Z; Z2];
    
    [X,IA] = unique(X,'rows');
    Z = Z(IA);
    
    %% Stats
    if mod(t,modtick) == 0
        s = 0;
        for i=1:size(X,1)
            s = s+Z{i}.f;
        end
        stats.mean_score(stati) = s/size(X,1);
    end
end
Ze = zeros(size(X,1),1);
Z2 = zeros(size(X,1),1);

for i=1:size(X,1)
    Ze(i) = Z{i}.r/Z{i}.n;
    Z2(i) = opt.fitness(X(i,:));
end
Z=Z2;
stats.mean_score(end) = mean(Z);
stats.mean_reward = stats.rewards ./ TN;
stats.expected_rewards = stats.expected_rewards ./ TN;
end

function z = stat_struct(f)
z = struct('r',0,'n',0,'f',f,'phat',0,'ci',[0 1]);
end

function v = evaluation(opt,type,z,Z)
p = sum(([Z.r]+1)./([Z.n]+2)) / length(Z);
n = length(Z);
ap = n*p;
bp = n-ap;
if strcmp(type,'greedy')
    v = z.r/z.n;
elseif strcmp(type,'thompson') || strcmp(type,'thompson_optimistic')
    exp = opt.exploration;
    a = z.r+ap;
    b = (z.n-z.r)+bp;
    v = betarnd(a/exp,b/exp);
    if strcmp(type,'thompson_optimistic')
        v = max(v,(z.r+ap)/(z.n+(ap+bp)));
    end
elseif strcmp(type,'epsilon_greedy')
    if rand() < opt.epsilon
        v = rand();
    else
        v = z.r/z.n;
    end
elseif strcmp(type,'ucb')
    v = (z.r+ap)/(z.n+ap+bp) + opt.exploration*sqrt(2*log(sum([Z.n])/(z.n+ap+bp)));
else
    error(strcat('unknown evaluation ', type))
end
end
