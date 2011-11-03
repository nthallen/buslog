function check_dir( route )
% check_dir(route)
r = almanac('earth','radius','sm');
for i = 1:length(route.direction)
    lat = route.direction(i).stops(:,2);
    lon = route.direction(i).stops(:,3);
    dist = distance(lat(1:end-1),lon(1:end-1),lat(2:end),lon(2:end)) ...
        * (pi/180) * r;
    fprintf(1,'Total distance for direction %d is %f miles\n',i,sum(dist));
    plot( lon, lat );
    hold on;
end
% for i = 1:length(route.path)
%     p = route.path(i).points;
%     plot(p(:,2),p(:,1));
% end
hold off;
