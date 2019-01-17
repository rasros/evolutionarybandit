%% Load
X_4C_8x10
opt = layoutopt(targets,weights);

%% Simulation

batchR = [];
itrR = [];
batchD = [];
batchRmax = [];
itrRmax = [];

%T is based on this calculation:
power_calculation(0.05,0.7,0.15,0.02)
opt.T = 4000;

opt.population = 20;
opt.mutation = 0.2;
opt.elitism = 8;

opt.exploration = 8;
opt.creation = 0.01;
opt.alpha = 0.3;

opt.repetitions = 20;
opt.crossover_rate = 0.9;

for k=1:1000
    ix = randperm(size(X0,1));
    X1 = X0(ix(1:opt.population),:);
    fs = [];
    for i=1:size(X1,1)
        fs = [fs opt.fitness(X1(i,:))];
    end
    [~,~,~,itrT,stats] = banditgasimulation(opt,X1);
    itrT = [0 itrT];
    itrR = [itrR; [mean(fs)/size(X1,1) stats.expected_rewards]];
    itrRmax = [itrRmax; [max(fs) stats.max_score]];

    ix = randperm(size(X0,1));
    X1 = X0(ix(1:opt.population),:);
    fs = [];
    for i=1:size(X1,1)
        fs = [fs opt.fitness(X1(i,:))];
    end
    [~,~,~,batchT,stats] = batchgasimulation(opt,X1);
    batchT = [0 batchT];
    batchR = [batchR; [mean(fs)/size(X1,1) stats.expected_rewards]];
    batchD = [batchD; [0 stats.duplicates]];
    batchRmax = [batchRmax; [max(fs) stats.max_score]];
    k
end

%% Save
bandit_results = struct();
bandit_results.R = itrR';
bandit_results.T = itrT';
bandit_results.max = itrRmax';

batch_results = struct();
batch_results.R = batchR';
batch_results.T = batchT';
batch_results.duplicates = batchD';
batch_results.max = batchRmax';

writetable(struct2table(bandit_results),'bandit_simulation_results2.csv')
writetable(struct2table(batch_results),'batch_simulation_results2.csv')

%% Calculate a-test
median(sum(itrR.*200,2))
iqr(sum(itrR.*200,2))
median(sum(batchR.*400,2))
iqr(sum(batchR.*400,2))

median(batchR(:,end))
iqr(batchR(:,end))
median(itrR(:,end))
iqr(itrR(:,end))

median(batchRmax(:,end))
iqr(batchRmax(:,end))
median(itrRmax(:,end))
iqr(itrRmax(:,end))

a_test(itrRmax(:,end), batchRmax(:,end))
a_test(itrR(:,end), batchR(:,end))
a_test(sum(itrR.*200,2), sum(batchR.*400,2))

%% Plots
%plot traces using same number of data points
f1=figure;
plot(batch_results.T,bandit_results.R(1:2:end,1:50), '-r');
ylim([0 0.4])
xlabel('Time [t]')
ylabel('Fitness [P(reward)]')

f2=figure;
plot(batch_results.T,batch_results.R(:,1:50), '-r');
ylim([0 0.4])
xlabel('Time [t]')
ylabel('Fitness [P(reward)]')

f3=figure;
plot(bandit_results.T,smooth(mean(bandit_results.R,2),5,'sgolay',3));
hold on
plot(batch_results.T,smooth(mean(batch_results.R,2),5,'sgolay',3),'--');
legend('Bandit', 'Baseline','Location','NorthWest')
xlim([0 opt.T]);
ylim([0 0.4]);
xlabel('Time [t]')
ylabel('Fitness [P(reward)]')

f4=figure;
cumbanditregret = bandit_results.T-mean(bandit_results.R,2).*bandit_results.T;
cumbatchregret = batch_results.T-mean(batch_results.R,2).*batch_results.T;
semilogx(bandit_results.T,cumbanditregret.*0.0001)
hold on
semilogx(batch_results.T,cumbatchregret.*0.0001,'--')
ylim([0 0.4])
xlim([0 4000])
xlabel('Time [t]')
ylabel('Cumulative regret [10 k]')
legend('Bandit', 'Baseline','Location','NorthWest')

%% matlab 2 tikz

matlab2tikz('f1.tex', 'figurehandle',f1,'width', '\fwidth', 'height', '\fheight')
matlab2tikz('f2.tex', 'figurehandle',f2,'width', '\fwidth', 'height', '\fheight')
matlab2tikz('f3.tex', 'figurehandle',f3,'width', '\fwidth', 'height', '\fheight')
matlab2tikz('f4.tex', 'figurehandle',f4,'width', '\fwidth', 'height', '\fheight')
