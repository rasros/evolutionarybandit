target = [3 2 2 2 2 2 2 2 2;
    3 3 4 4 1 1 1 4 4;
    3 3 4 4 1 1 1 4 4;
    3 3 1 1 1 1 1 1 1;
    3 3 1 1 1 1 1 1 1;
    3 3 1 1 4 4 4 1 1;
    3 1 1 1 4 4 4 1 1;];
target = flip(target);

%% Plot grid
%Cmap = [254,254,246;
%    202, 1, 7;
%    52, 73, 128;
%    255 215, 3]./255;
Cmap = [254,254,246,153
        102,194,165,153;
        252,141,98,153;
        141,160,203,153;]./255;
 
figure
hold on
axis equal
axis off
grid off
%xlim([-1 10])
%ylim([-1 8])
for i=1:size(target,1)
    for j=1:size(target,2)
        rectangle('Position', [j,i,1,1], ...
            'FaceColor', Cmap(target(i,j),:), ...
            'EdgeColor', 'k');
    end
end
matlab2tikz('dissection1.tex', 'width', '\fwidth', 'height', '\fheight')

%% Plot finished result
R = grid2rect(target);
plotrect(R,'paper',4,1.5);
grid off
%xlim([-1 10])
%ylim([-1 8])
matlab2tikz('dissection3.tex', 'width', '\fwidth', 'height', '\fheight')

%% Plot intersetcions

C1 = [2 1;
     1 1;
     1 0;
     4 0;
     4 2;
     7 2;
     7 0;
     9 0;
     9 4;
     7 4;
     7 6;
     4 6;
     4 4;
     2 4;];
 
 C2 = [0 0;
       0 7;
       1 7;
       1 6;
       2 6;
       2 1;
       1 1;
       1 0];

figure
hold on
axis equal
axis off
grid off
rectangle('Position', [4,0,3,2], ...
    'FaceColor', Cmap(4,:), ...
    'EdgeColor', 'k');
rectangle('Position', [7,4,2,2], ...
    'FaceColor', Cmap(4,:), ...
    'EdgeColor', 'k');
rectangle('Position', [2,4,2,2], ...
    'FaceColor', Cmap(4,:), ...
    'EdgeColor', 'k');
rectangle('Position', [1,6,8,1], ...
    'FaceColor', Cmap(2,:), ...
    'EdgeColor', 'k');

fill(C1(:,1),C1(:,2),Cmap(1,1:3),'FaceAlpha', 0.6);
fill(C2(:,1),C2(:,2),Cmap(3,1:3),'FaceAlpha', 0.6);
line([4 4], [2 4], 'LineStyle', ':', 'Color', 'k', 'LineWidth', 2);
line([7 7], [2 4], 'LineStyle', ':', 'Color', 'k');
line([4 7], [4 4], 'LineStyle', ':', 'Color', 'k');
line([1 1], [1 6], 'LineStyle', ':', 'Color', 'k');
plot(4,2,'bo','MarkerSize',6)
plot(2,1,'bo','MarkerSize',6)
plot(4,4,'bo','MarkerSize',6)
plot(7,2,'bo','MarkerSize',6)
plot(7,4,'bo','MarkerSize',6)
plot(1,1,'ro','MarkerSize',6)
plot(1,6,'ro','MarkerSize',6)
%xlim([-1 10])
%ylim([-1 8])


matlab2tikz('dissection2.tex', 'width', '\fwidth', 'height', '\fheight')
