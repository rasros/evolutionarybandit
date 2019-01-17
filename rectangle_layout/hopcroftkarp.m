function [P] = hopcroftkarp(U,E)
%HOPCROFTKARP Algorithm adapted from wikipedia pseudo code.
%    Return value P is the maximum matching.

N = size(E,1);
NIL = N+1;
D = zeros(N+1,1);

Pair_U = ones(N,1)*NIL;
Pair_V = ones(N,1)*NIL;

    function res = BFS()
        Q = [];
        for uui = 1:length(U)
            uu = U(uui);
            if Pair_U(uu) == NIL
                D(uu) = 0;
                Q = [Q uu];
            else
                D(uu) = Inf;
            end
        end
        D(NIL) = Inf;
        while ~isempty(Q)
            u = Q(1);
            Q(1) = [];
            if D(u) < D(NIL)
                %Adj = find(E(u,:)==1);
                Adj = E{u};
                for vi = 1:length(Adj)
                    v = Adj(vi);
                    if D(Pair_V(v)) == Inf
                        D(Pair_V(v)) = D(u) + 1;
                        Q = [Q Pair_V(v)];
                    end
                end
            end
        end
        res = D(NIL) ~= Inf;
    end

    function res = DFS(u)
        if u ~= NIL
            %Adj = find(E(u,:) == 1);
            Adj = E{u};
            for vi = 1:length(Adj)
                v = Adj(vi);
                if D(Pair_V(v)) == D(u) + 1
                    if DFS(Pair_V(v))
                        Pair_V(v) = u;
                        Pair_U(u) = v;
                        res = 1;
                        return;
                    end
                end
            end
            D(u) = Inf;
            res = 0;
            return;
        end
        res = 1;
    end

M = [];
while BFS()
    for ui=1:length(U)
        if Pair_U(U(ui)) == NIL
            if DFS(U(ui))
                M = [M U(ui)];
            end
        end
    end
end

P = [M; Pair_U(M)'];
end