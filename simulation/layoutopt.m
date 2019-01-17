function opt = layoutopt(targets,weights)

%% Default parameters
opt = struct();
[n,m,d] = size(targets);

C = max(targets(:));
opt.C = C;
R = zeros(d,1);
for i=1:d
    R(i) = size(grid2rect(targets(:,:,i),'graph',C),1);
end
opt.T = 4000;
opt.N = n*m;

%% Hyper parameters

opt.population = 5;
opt.elitism = 2;

opt.mutation = 0.5;

%only used by batch
opt.repetetions = 10;
opt.crossover_rate = 0.8;

%only used by iterative
opt.alpha = 0.4;
%opt.evaluation = 'greedy';
%opt.evaluation = 'thompson';
%opt.evaluation = 'thompson_optimistic';
opt.evaluation = 'ucb';
opt.exploration = 10;
opt.creation = 0.1;

%opt.crossover_type = 'scatter';
%opt.crossover_type = 'single_point';
opt.crossover_type = 'rectangle';

opt.experiment = @(f,n) binornd(n,f);
opt.generate = @() generate(C,n,m);
opt.fitness = @(X) fitness(n,m,X,targets,weights,C,R);
opt.mutate = @(G) mutate(G,C,n,m);
opt.crossover = @(G1,G2) crossover(opt,G1,G2,n,m);
opt.regenerate = @(X) regenerate(opt,X);
end

%% Genetic operators
function f = fitness(n,m,G,targets,weights,C,R)
x = reshape(G,n,m);
fs = ones(1,length(weights));
for i=1:length(weights)
    neq = sum(sum(double(targets(:,:,i)~=x)));
    reg = abs(R(i) - size(grid2rect(x,'greedy',C),1));
    fs(i) = 1-min(n*m,neq+reg)/(n*m);
    fs(i) = weights(i)*fs(i);
end
f = 1-harmmean(1-fs);
end

function G = generate(C,n,m)
G = gridrnd(C,n,m);
G = reshape(G,1,n*m);
end

function x = regenerate(opt,X)
if strcmp(opt.regeneration, 'mutate')
    ix = X(randi(size(X,1), 1));
    x = opt.mutate(X(ix,:));
elseif strcmp(opt.regeneration, 'crossover')
    ix1 = randi(size(X,1));
    ix2 = ix1;
    while ix1==ix2
        ix2 = randi(size(X,1));
    end
    x = opt.crossover(X(ix1,:), X(ix2,:));
elseif strcmp(opt.regeneration, 'never')
    x = X(randi(size(X,1)),:);
elseif strcmp(opt.regeneration, 'generate')
    x = opt.generate();
else
    error('unknown opt.regenerate')
end
end

function G = mutate(G,C,n,m)
G = reshape(G,n,m);
R = grid2rect(G,'greedy',C);
ix = randi(size(R,1));
r = R(ix,:);
Cp = setdiff(1:C,r(1));
c = Cp(randi(C-1));
x1 = r(2)+randi(r(4));
x2 = x1+randi(r(2)+r(4)-x1+1)-1;
y1 = r(3)+randi(r(5));
y2 = y1+randi(r(3)+r(5)-y1+1)-1;
G(y1:y2,x1:x2) = c;
G = reshape(G,1,n*m);
end

function [GS1,GS2] = crossover(opt,G1,G2,n,m)
%CROSSOVER Produce crossovers between G1 and G2.
if strcmp(opt.crossover_type, 'scatter')
    GS1 = zeros(1,n*m);
    GS2 = zeros(1,n*m);
    for i=1:n*m
        if rand() > 0.5
            GS1(i) = G1(i);
            GS2(i) = G2(i);
        else
            GS1(i) = G2(i);
            GS2(i) = G1(i);
        end
    end
elseif strcmp(opt.crossover_type, 'single_point')
    G1 = reshape(G1,n,m);
    G2 = reshape(G2,n,m);
    if rand() > 0.5
        x = randi(m-1);
        GS1 = [G1(:,1:x) G2(:,x+1:end)];
        GS2 = [G2(:,1:x) G1(:,x+1:end)];
    else
        y = randi(n-1);
        GS1 = [G1(1:y,:); G2(y+1:end,:)];
        GS2 = [G2(1:y,:); G1(y+1:end,:)];
    end
    if rand() > 0.5
        tmp = GS1;
        GS1 = GS2;
        GS2 = tmp;
    end
    GS1 = reshape(GS1,1,n*m);
    GS2 = reshape(GS2,1,n*m);
elseif strcmp(opt.crossover_type, 'rectangle')
    G1 = reshape(G1,n,m);
    G2 = reshape(G2,n,m);
    w = m+1;
    h = n+1;
    while w > m || w <= 0  
        w = round(normrnd(sqrt(m), 2*m));
    end
    while h > n || h <= 0  
        h = round(normrnd(sqrt(n), 2*n));
    end
    x = randi(m-w+1);
    y = randi(n-h+1);
    mask = false(n,m);
    mask(y:y+h-1,x:x+w-1) = true;
    GS1 = G1.*(G1&mask) + G2.*(G2&~mask);
    GS2 = G1.*(G1&~mask) + G2.*(G2&mask);
    if rand() > 0.5
        tmp = GS1;
        GS1 = GS2;
        GS2 = tmp;
    end
    GS1 = reshape(GS1,1,n*m);
    GS2 = reshape(GS2,1,n*m);
end
end