function route = route_def(route_num)
% route = route_def(route_num);
if route_num == 77
    route = route77;
elseif route_num == 78
    route = route78;
else
    error('Unknown route number');
end