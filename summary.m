function summary(pattern, dirs, xl)
% summary(pattern, dir, xl);
% This is the graph I've been using (1/15)
% Note that plotByDirection has a hard-coded dst value that must be
% adjusted.
if nargin < 2
    dirs = [];
end
if nargin < 3
    xl = [];
end
files = dir(pattern);
for i=1:length(files)
    A = work_path(files(i).name);
    ax = plotByDirection(A,dirs);
    if ~isempty(xl)
        set(ax,'xlim',xl);
    end
end