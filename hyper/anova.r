data = read.csv("random_search_bandit.csv")
attach(data)
av1 = aov(R~population*mutation*elitism*alpha*exploration*creation)
av2 = aov(R~population*mutation*elitism*repetitions*crossover_rate)
summary(av1)
summary(av2)
detach(data)
