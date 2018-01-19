function plotSpeedByDistance(A)
rundirs = unique([A.run A.newdirs],'rows');
for dir = 1:length(A.route.direction)
    figure;
    v = find(rundirs(:,2)==dir);
    for run = rundirs(v,1)'
        x = find(A.run == run);
        dd = diff(A.dist(x));
        dt = diff(A.time(x));
        xx = find(dt > .002);
        plot(A.dist(x(xx)),dd(xx)./dt(xx),'.-');
        hold on;
    end
    hold off;
    title([A.route.route ': ' A.route.direction(dir).title ' ' ...
        A.route.direction(dir).name ]);
end