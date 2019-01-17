function [R,t1,t2] = dissectpoly(G)
%DISSECTPOLY dissects orthogonal polygon with holes. Grid G is a binary
%matrix where 1 represent parts of the polygon. There can be multiple
%polygons. The dissection works by finding diagonals that connect two 
%concave corners that are axis aligned and forming a bipartite graph of the
%horizontal and vertical diagonals. The edges in the graph represent
%intersections of the diagonals. The maximum independent set is then used
%to find diagonals that are used in the dissection. We can find maximum
%independent set using Hopcroft-Karp to obtain maximum matching and then
%using Konigs theorem to generate maximum independent set. This function is
%fairly optimized and can solve matrices of eg. 1E4x1E4 in seconds.

%#ok<*AGROW>
 
%Create bipartite graph by connecting nodes with axis aligned diagonals
tic;
[n,m] = size(G);
%Padding grid for easier indexing.
G = [zeros(n+2,1) [zeros(1,m); G; zeros(1,m)] zeros(n+2,1)];
[n,m] = size(G);


%% Forming graph of intersecting cuts

%Find concave corners, which is any 2x2 region with one 0 and three 1.
%It is further divided into corners that are oriented approprietly.
%Only corners with south east or south west rotation can be a starting
%point for a vertical line, and south east and north east for
%horizontal.

Cn = 0;
CI = zeros(16,1,'uint32');
CJ = zeros(16,1,'uint32');
C1I = zeros(16,1,'uint32');
C2I = zeros(16,1,'uint32');
C3I = zeros(16,1,'uint32');
C4I = zeros(16,1,'uint32');
C1J = zeros(16,1,'uint32');
C2J = zeros(16,1,'uint32');
C3J = zeros(16,1,'uint32');
C4J = zeros(16,1,'uint32');
C1IX = zeros(16,1,'uint32');
C2IX = zeros(16,1,'uint32');
C3IX = zeros(16,1,'uint32');
C4IX = zeros(16,1,'uint32');
C1n = 0;
C2n = 0;
C3n = 0;
C4n = 0;
for j=2:m-1
    for i=2:n-1
        if G(i,j)+G(i+1,j)+G(i,j+1)+G(i+1,j+1) == 3
            Cn = Cn+1;
            CI(Cn) = i;
            CJ(Cn) = j;
            if Cn == size(CI,1)+1
                CI = [CI; zeros(size(CI,1)*2,1,'uint32')];
                CJ = [CJ; zeros(size(CJ,1)*2,1,'uint32')];
            end
            if G(i,j) == 0
                C1n = C1n+1;
                if C1n == size(C1I,1)+1
                    C1I = [C1I; zeros(size(C1I,1)*2,1,'uint32')];
                    C1J = [C1J; zeros(size(C1J,1)*2,1,'uint32')];
                    C1IX = [C1IX; zeros(size(C1IX,1)*2,1,'uint32')];
                end
                C1I(C1n) = i;
                C1J(C1n) = j;
                C1IX(C1n) = Cn;
            end
            if G(i,j+1) == 0
                C2n = C2n+1;
                if C2n == size(C1I,1)+1
                    C2I = [C2I; zeros(size(C2I,1)*2,1,'uint32')];
                    C2J = [C2J; zeros(size(C2J,1)*2,1,'uint32')];
                    C2IX = [C2IX; zeros(size(C2IX,1)*2,1,'uint32')];
                end
                C2I(C2n) = i;
                C2J(C2n) = j;
                C2IX(C2n) = Cn;
            end
            if G(i+1,j) == 0
                C3n = C3n+1;
                if C3n == size(C1I,1)+1
                    C3I = [C3I; zeros(size(C3I,1)*2,1,'uint32')];
                    C3J = [C3J; zeros(size(C3J,1)*2,1,'uint32')];
                    C3IX = [C3IX; zeros(size(C3IX,1)*2,1,'uint32')];
                end
                C3I(C3n) = i;
                C3J(C3n) = j;
                C3IX(C3n) = Cn;
            end
            if G(i+1,j+1) == 0
                C4n = C4n+1;
                if C4n == size(C1I,1)+1
                    C4I = [C4I; zeros(size(C4I,1)*2,1,'uint32')];
                    C4J = [C4J; zeros(size(C4J,1)*2,1,'uint32')];
                    C4IX = [C4IX; zeros(size(C4IX,1)*2,1,'uint32')];
                end
                C4I(C4n) = i;
                C4J(C4n) = j;
                C4IX(C4n) = Cn;
            end
        end
    end
end
CI(Cn+1:end) = [];
CJ(Cn+1:end) = [];
C1I(C1n+1:end) = [];
C1J(C1n+1:end) = [];
C1IX(C1n+1:end) = [];
C2I(C2n+1:end) = [];
C2J(C2n+1:end) = [];
C2IX(C2n+1:end) = [];
C3I(C3n+1:end) = [];
C3J(C3n+1:end) = [];
C3IX(C3n+1:end) = [];
C4I(C4n+1:end) = [];
C4J(C4n+1:end) = [];
C4IX(C4n+1:end) = [];
C = [CJ CI];

% U and V are the disjoint parts of the bipartite graph, which are
% vertical and horizontal diagonals respectively.
U = [];
V = [];
nn = 0;

% Rows in L are the diagonals connecting two concave corners, i.e. cuts in
% the polygon.
L = [];

% Sorting these to allow faster construction of axis parallell lines
U1 = [C1J C1I C1IX; C2J C2I C2IX];
U2 = sortrows([C3J C3I C3IX; C4J C4I C4IX], [1,2]);
V1 = [C1J C1I C1IX; C3J C3I C3IX];
V2 = sortrows([C2J C2I C2IX; C4J C4I C4IX], [2,1]);

%Vertical lines, from 1/2 type corners to 3/4
for i=1:(C1n+C2n)
    n1 = U1(i,:);
    [x1,x2] = findinsorted(U2(:,1), n1(1));
    for j=1:x2-x1+1
        n2 = U2(x1+j-1,:);
        if n1(2) < n2(2)
            free = 1;
            for ei=n1(2)+1:n2(2)
                if G(ei,n1(1)) == 0 || G(ei,n1(1)+1) == 0
                    free = 0;
                    break
                end
            end
            if free
                nn = nn+1;
                U = [U; nn];
                L(nn,:) = [n1(3); n2(3)];
            end
        end
    end
end
%Horizontal lines, from 1/3 type corners to 2/4
for i=1:(C1n+C3n)
    n1 = V1(i,:);
    [y1,y2] = findinsorted(V2(:,2),n1(2));
    for j=1:y2-y1+1
        n2 = V2(y1+j-1,:);
        if n1(1) < n2(1)
            free = 1;
            for ei=n1(1)+1:n2(1)
                if G(n1(2),ei) == 0 || G(n1(2)+1,ei) == 0
                    free = 0;
                    break
                end
            end
            if free
                nn = nn+1;
                V = [V; nn];
                L(nn,:) = [n1(3); n2(3)];
            end
        end
    end
end

% E is the edges in the graph.
E = cell(nn,1);

% Find intersection, which become edges in the graph

if ~isempty(V) && ~isempty(U)
    ys = sortrows([C(L(V,1),2) V],1);
    for i=1:length(U)
        u = U(i,:);
        x = C(L(u,1),1);
        y1 = C(L(u,1),2);
        y2 = C(L(u,2),2);
        [yr1,yr2] = findinsorted(ys(:,1),y1:y2);
        for j=yr1:yr2
            v = ys(j,2);
            x1 = C(L(v,1),1);
            x2 = C(L(v,2),1);
            if x1 <= x && x2 >= x
                E{u} = [E{u}; v];
                E{v} = [E{v}; u];
            end
        end
    end
end
t1 = toc;
tic;

%% Calculating included cuts
% Find maximum matching P using Hopcroft-Karp, with bipartite graph parts U
% and V, and edges defined in E.

P = hopcroftkarp(U,E);

% Convert to minimum vertex cover and then to maximum independent set, by
% following an alternating path starting from unmatched vertices, i.e.
% Konigs algorithm.
mincover = zeros(1,nn);
    function visit_unmatched(node)
        mincover(node) = -1;
        adj = E{node};
        rec = [];
        for ia=1:length(adj)
            if mincover(adj(ia)) == 0
                mincover(adj(ia)) = 1;
                [ii,jj] = find(P == adj(ia));
                mate = P(2-(ii-1),jj);
                rec = [rec; mate];
            end
        end
        for ir=1:length(rec)
            visit_unmatched(rec(ir));
        end
    end

if isempty(P)
    unmatched = U;
else
    unmatched = setdiff(U,P(1,:));
end
for i=1:length(unmatched)
    visit_unmatched(unmatched(i))
end
Z = find(mincover ~= 0);
if isempty(Z)
    maxindset = [];
else
    K = union(setdiff(U,Z), intersect(V,Z));
    %Maximum independent set is the edges we must use in the dissection.
    maxindset = setdiff(1:nn,K);
end
t2 = toc;
tic;

%% Convert cuts to rectangles
% We create a horizontal and vertical mask which is passed on to
% gridcuts2rect.
horizmask = false(n,m);
vertmask = false(n,m);
for i=1:length(maxindset)
    e = maxindset(i);
    n1 = C(L(e,1),:);
    n2 = C(L(e,2),:);
    if n1(1) == n2(1)
        vertmask(n1(2)+1:n2(2), n1(1)+1:n2(1)+1) = true;
    else
        horizmask(n1(2)+1:n2(2)+1, n1(1)+1:n2(1)) = true;
    end
end
R = gridcuts2rect(G,horizmask,vertmask);

end
