function axout = plotByDirection(A, dirs)
% ax = plotByDirection(A[, dirs]);
if nargin < 2 || isempty(dirs)
    dirs = 1:length(A.route.direction);
end
dst = 0;
colors = 'bgry';
fig = figure;
hda = datacursormode(fig);
set(hda, 'UpdateFcn', @buslog_cursor_text_func);
datacursormode(fig, 'off');
ax = gca;
legends = cell(length(dirs),1);
handles = zeros(length(dirs),1);
for dirnum = 1:length(dirs)
    dir = dirs(dirnum);
    x = A.route.departures{dir};
    y = zeros(size(A.route.departures{dir}));
    handles(dirnum) = ...
        plot(x, y, [ 'o' colors(dirnum) ]);
    hold on;
    legends{dirnum} = sprintf('%s %s', A.route.direction(dir).title, ...
        A.route.direction(dir).name );
end
for dirnum = 1:length(dirs)
    dir = dirs(dirnum);
    for run = find([A.runs.dir] == dir)
        plot(A.runs(run).time+dst,A.runs(run).dist,['.-' colors(dirnum)]);
    end
end
hold off;
title([A.route.route ': ' A.day ]);
legend( handles, legends);
if nargout > 0
    axout = ax;
end
