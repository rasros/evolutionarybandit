function plotrect(R,map,C,linewidth)
figure
hold on
axis equal
axis off
grid on
view(0,-90)
if nargin <= 1
    map = '';
end
if nargin <= 2
    C = max(R(:,1));
end
if nargin <= 3
    linewidth = 3;
end
if strcmp(map,'piet')
    Cmap = [254,254,246;
         202, 1, 7;
         52, 73, 128;
         255 215, 3]./255;
     C = 4;
elseif strcmp(map, 'paper')
    Cmap = [254,254,246,153
        102,194,165,153;
        252,141,98,153;
        141,160,203,153;]./255;
else
    Cmap = parula(double(C));
end
for i = 1:size(R,1)
    c = mod(R(i,1),C);
    if c == 0
        c = C;
    end
    x = R(i,2);
    y = R(i,3);
    w = R(i,4);
    h = R(i,5);
    c = Cmap(c,:);
    rectangle('Position', [x,y,w,h], ...
        'FaceColor', c, ...
        'EdgeColor', 'k',...
        'LineWidth', linewidth);
end
end