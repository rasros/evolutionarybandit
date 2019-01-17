function G = gridrnd(C,n,m)
%GRIDRND Generates random grid.
%   C: number of colors
%   n,m: rows and columns

G = zeros(n,m);
while find(G==0,1)
    ix = find(G==0);
    ix = ix(randi(length(ix)));
    [y1,x1] = ind2sub([n,m],ix);
    y2 = y1+randi(n-y1+1)-1;
    x2 = x1+randi(m-x1+1)-1;
    if find(G(y1:y2,x1:x2) > 0)
        %TODO this naive algo can be improved
        Gs = G(y1:y2,x1:x2);
        bestarea = 1;
        bestj = 1;
        besti = 1;
        for i=1:size(Gs,1)
            for j=1:size(Gs,2)
                if i*j <= bestarea
                    continue
                end
                if isempty(find(Gs(1:i,1:j)>0,1))
                    bestj = j;
                    besti = i;
                    bestarea = i*j;
                end
            end
        end
        x2 = x1+bestj-1;
        y2 = y1+besti-1;
    end
    G(y1:y2,x1:x2) = randi(C);
end
