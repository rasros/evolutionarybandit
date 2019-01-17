function R = grid2rect(G,method,C)
if nargin == 1
    method = 'graph';
end
if nargin <= 2
    C = double(max(G(:)));
end
R = uint32([]);
[n,m] = size(G);
for c=1:C
    Gc = G == c;
    if strcmp(method, 'jacop')
        used = logical(1-Gc);
        for i=1:n
            for j=1:m
                if used(i,j) == 0
                    Rr = region(Gc,j,i);
                    minx = min(Rr(:,1));
                    miny = min(Rr(:,2));
                    maxx = max(Rr(:,1));
                    maxy = max(Rr(:,2));
                    Gcs = accumarray(Rr(:,1:2), 1)';
                    Gcs = Gcs(miny:maxy,minx:maxx);
                    used(miny:maxy,minx:maxx) = used(miny:maxy,minx:maxx) | Gcs;
                    [ns,ms] = size(Gcs);
                    if ns*ms == sum(Gcs(:))
                        Rc = [c minx-1 miny-1 maxx-minx+1 maxy-miny+1];
                    else
                        Rcs = piet.Rectangles.rectangles(Gcs);
                        Rc = [ones(size(Rcs,1),1)*c Rcs(:,1)+minx-1 Rcs(:,2)+miny-1 Rcs(:,3:4)];
                    end
                    R = [R; Rc];
                end
            end
        end 
    elseif strcmp(method, 'graph')
        Rc = dissectpoly(Gc);
        R = [R; ones(size(Rc,1),1,'uint32')*c Rc];
    elseif strcmp(method, 'greedy')
        Gc0 = [zeros(n+2,1) [zeros(1,m); Gc; zeros(1,m)] zeros(n+2,1)];
        Rc = gridcuts2rect(Gc0,zeros(n+2,m+2),zeros(n+2,m+2));
        R = [R; ones(size(Rc,1),1,'uint32')*c Rc];
    else
         error('bad method type')
    end
end
end

function [R,free] = region(G,x,y,free)
% All contiguous regions from initial coordinates
[n,m] = size(G);
if nargin == 3
    free = true([n,m]);
end
assert(G(y,x))
R = int8([x,y]);
free(y,x) = 0;
if x<m && G(y,x+1) == 1 && free(y,x+1)
    [R2,free] = region(G,x+1,y,free);
    R = [R; R2];
end
if x-1>0 && G(y,x-1) == 1 && free(y,x-1)
    [R2,free] = region(G,x-1,y,free);
    R = [R; R2];
end
if y<n && G(y+1,x) == 1 && free(y+1,x)
    [R2,free] = region(G,x,y+1,free);
    R = [R; R2];
end
if y-1>0 && G(y-1,x) == 1 && free(y-1,x)
    [R2,free] = region(G,x,y-1,free);
    R = [R; R2];
end
end

