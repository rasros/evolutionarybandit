function R = gridcuts2rect(G,hG,vG)
%GRID_CUTS_TO_RECTANGLES Creates rectangles from a grid cut horizontally
%and vertically. The cuts are defined as a mask which are stopping points
%when iterating through the graph. Expects zero padded grid.
%   G: grid which defines a polygon, 1s are part of the polygon and 0s are
%   outside.
%   hG: horizontal cuts in the polygon.
%   vG: vertical cuts in the polygon.
%#ok<*AGROW>

[n,m] = size(G);
vG = vG | ~G;

nr = 0;
X = zeros(256,1,'uint32');
Y = zeros(256,1,'uint32');
W = zeros(256,1,'uint32');
H = zeros(256,1,'uint32');

[I,J] = find(G);

for ix=1:length(I)
    i=I(ix);
    j=J(ix);
    if ~G(i,j)
        % Need to check this dynamically since G is changed during
        % iterations.
        continue
    end
    
    % grow vertically with width 1
    for h = 1:n-i
        if ~G(i+h,j) || hG(i+h,j) > 0
            break
        end
        G(i+h,j) = 0;
    end
    
    % grow horizontally with height h, we don't need to check G here
    % since it is included in vertmask, and there cannot be any
    % rectangles to the right.
    for w = 1:m-j
        walled = 0;
        for wi=i:i+h-1
            if vG(wi,j+w) > 0
                walled = 1;
                break
            end
        end
        if walled
            break
        end
        G(i:i+h-1,j+w) = 0;
    end

    nr = nr+1;
    if nr == size(X,1)+1
        X = [X; zeros(size(X,1)*2,1,'uint32')];
        Y = [Y; zeros(size(X,1)*2,1,'uint32')];
        W = [W; zeros(size(X,1)*2,1,'uint32')];
        H = [H; zeros(size(X,1)*2,1,'uint32')];
    end
    X(nr) = j-2;
    Y(nr) = i-2;
    W(nr) = w;
    H(nr) = h;
end
R = [X(1:nr) Y(1:nr) W(1:nr) H(1:nr)];
end
