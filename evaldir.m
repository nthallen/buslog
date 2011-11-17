function [ dist, xdist ] = evaldir( A, route, ix )
% seg = evaldir( A, route, ix, dir );
% A is the bus location structure
% route is the route definition structure
% ix is a vector of indices into the bus location structure
% dir is the index of the direction against which to evaluate
%
% Returns 3 matrices of size (length(ix),length(route.direction))
% dist distance between readings along route
% xdist distance from the route
np = length(ix);
nd = length(route.direction);
dist = zeros(np-1, nd);
xdist = zeros(np, nd);
rfact = almanac('earth','radius','sm')*2*pi/360; % miles per degree
for dir = 1:length(route.direction)
    P = route.direction(dir).stops(:,[2 3]);
    rdist = route.direction(dir).dist;
    rcdist = route.direction(dir).cdist;
    for j = 1:length(ix)
        rxdist = distance([A.lat(ix(j)) A.lon(ix(j))], P);
        st = find(rxdist == min(rxdist));
        % I'll assume only one.
        st = st(1);
        % Given P1, P2 on route and bus position B
        % compare dist(P1,B)+dist(B,P2) with dist(P1,P2)
        % (d1+d2-d0)/d0 and choose the one closer to zero
        if st ~= 1 && st == length(rxdist)
            st = st - 1;
        elseif st ~= 1
            q0 = rxdist(st)+rxdist(st-1)-rdist(st-1);
            q1 = rxdist(st)+rxdist(st+1)-rdist(st);
            if q0 < q1
                st = st-1;
            end
        end
        dist(j,dir) = rcdist(st) + rxdist(st);
        xdist(j,dir) = (rxdist(st)+rxdist(st+1)-rdist(st))/2;
    end
end
dist = dist * rfact;
xdist = xdist * rfact;