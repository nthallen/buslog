%function work_path
% First, get the appropriate path and calculate distances
% Next, get bus locations. For each location:
%  determine which path segment the location is closest to
%  and the distance along that segment.
%
% For each direction in the route:
% Given path P(i), 1 <= i <= N
% for 1 <= j < N, calculate:
%   Angle between P(j) and P(j+1)
%   (distance from P(j) to P(j+1))
%   Cartesian coordinates for P(j)
%   Cartesian coordinates of normal vector N(j)
%   Cartesian coordinates of coplanar vector D(j) = P(j)xN(j)
%   Angle phi(j) = asin(dot(P(j+1),D(j)) (verify)
%
% load day's log
% for each vehicle, break data into directions
%   this is complicated when the dirTag is not a direction.
%   I will adopt the heuristic to replace unknown tags with
%   the most recent recognized tag. Or I could discard them
%   at the end or beginning of a run.
% for location C
%   for 1 <= j < N, calculate
%   Cartesian coordinates of coplanar vector E = ||CxN(j)||
%   Cartesian coordinates of coplanar vector C' = ||N(j)xE||
%   Angle zeta between C' and P(j) = asind(dot(C',D(j))
%   If 0 <= zeta <= phi(j), then we are adjacent to the
%   segment. Angle to the segment is acos(dot(C,N(j)))-pi/2
%   Otherwise, distance to the segment is distance to the
%   nearest endpoint.
route = route77;
% coslat = cosd(P(:,1));
% CP = [ coslat .* cosd(P(:,2))  coslat .* sind(P(:,2))  sind(P(:,1)) ];
% N = cross(CP(1:end-1,:),CP(2:end,:));
% N = diag(1./sqrt(sum(N.*N,2))) * N;
% D = cross(CP(1:end-1,:),N); % should be pre-normalized
% dist2 = abs(asind(dot(CP(2:end,:),D,2)));
% x = 1:size(N,1);
% plot(x,dist,'+',x,dist2,'*');
A = ReadBusLog('route_77/77_111107.log');
A.time = time2d(A.time)/3600-5;
rfact = almanac('earth','radius','sm')*2*pi/360; % miles per degree
%% Fill in direction information
dirs = zeros(size(A.time));
for dir = 1:length(route.direction)
    dirs(strcmp(route.direction(dir).tag,A.dirTag)) = dir;
    P = route.direction(dir).stops(:,[2 3]);
    route.direction(dir).dist = distance(P(1:end-1,:),P(2:end,:));
    route.direction(dir).cdist = [ 0; cumsum(route.direction(dir).dist) ];
end
%% Unknown directions have index of zero
buses = unique(A.ID);
for bus = 1:length(buses)
    isbus = find(strcmp(A.ID, buses{bus} ));
    % find unidentified directions. If they are preceeded
    % and followed by the same direction, replace them with
    % that direction.
    % unknowns are the starting index in isbus of each run of unknowns
    unknowns = find(diff(dirs(isbus) == 0) > 0);
    for unk = 1:length(unknowns)
        % kk(1) is the first index in isbus with a known direction
        kk = find(diff(dirs(isbus) ~= 0 & isbus > isbus(unknowns(unk))));
        if ~isempty(kk) && dirs(isbus(kk(1)+1)) == dirs(isbus(unknowns(unk)))
            % fill in
            dirs(isbus(unknowns(unk)+1:kk(1))) = dirs(isbus(kk(1)+1));
        end
    end
end
%% Parse bus route into runs and assign directions
% I think I need to look at a bus's entire route rather than one
% direction at a time. I should split a route into directions, but
% then analyze distance along direction to decide if the points are
% in the correct direction. Any backtracking is probably an argument
% for splitting. Then if there is a sequence of points that are
% backtracking, they probably belong on another direction. The obvious
% candidate direction would be the one that follows, but we could also
% try against other directions. For the routes I've looked at so far,
% the direction's 'name' attribute is Inbound or Outbound, so switching
% to the opposite of the current direction might be a good guess.
%
% 1 2 3 4 5 4 5 6 7 8 all Outbound could be split to:
% 1 2 3 4 5 Outbound
% 4 5 6 7 8 Outbound
%
% This scenario is perhaps more likely to be seen as:
%
% 1 2 3 4 1 2 3 4 Outbound to become
% 1 2 3 4 Outbound
% 1 2 3 4 Outbound
%
% 1 2 3 4 5 4 3 Outbound followed by
% 2 1 Inbound
%
% Split first into:
% 1 2 3 4 5 Outbound
% 5 4 3 Outbound but backwards
% 2 1 Inbound
%
% Then merge last two as Inbound
%
% Also possible:
% 1 2 3 4 5 4 3 2 1 2 3 4 5 Outbound =>
% 1 2 3 4 5 Outbound
% 5 4 3 2 1 Outbound but backwards
% 1 2 3 4 5 Outbound
%
% Try backwards section on other directions to find the best fit
%
% Besides going backwards, it is also possible for a bus to go off the
% route. This could happen due to diversion or due to simply reporting
% the wrong direction (Arlmont vs Heights). In the latter case, the
% entire run should be switched if it can be determined that another
% direction is a better fit.
%
% For comparing performance, the 78 has a common path for most of the run,
% so it would make sense to combine data for both inbound routes and both
% outbound routes where they overlap, and then provide separate graphs for
% the unique parts. That requires route analysis... But the route data
% gives a clue in that the stops on the different routes are identified as
% being the same for part of the run.
%
% I think the direction data from the buses is too unreliable. It probably
% makes sense to evaluate all points along all directions, and then choose
% the best fit, i.e. the longest run of matching points.
%
newdirs = dirs * 0;
n_runs = 0;
for bus = 1:length(buses)
    isbus = find(strcmp(A.ID, buses(bus) ));
    [ dist, xdist ] = evaldir(A, route, isbus);
    fwd = diff(dist) > -.01;
    tgap = (diff(A.time(isbus))>.2);
    fwd(tgap,:) = 0;
    % if the direction is right and the time gap is small enough,
    fwdfrom = [fwd; fwd(end,:)];
    fwd = [ fwd(1,:); fwd ];
    near = xdist < .02;
    good = fwd & near;
    goodfrom = fwdfrom & near;
    good = [ good; 0*good(end,:) ];
    run0 = 1;
    while run0 < length(isbus)
        % Find the next good point
        good(run0,:) = goodfrom(run0,:);
        run1 = min(find(any(good(run0:end,:)')));
        if isempty(run1)
            fprintf(1,'%s: Discarding %d points at end of run\n', ...
                buses{bus}, length(isbus)-run0+1);
            run0 = length(isbus);
        else
            if run1 > 1
                fprintf(1, '%s: Discarding %d points\n', buses{bus}, run1-1);
                run0 = run0 + run1 - 1;
            end
            % Now select the longest run
            gi = find(good(run0,:));
            bestdir = 0;
            bestlen = 1;
            for i = 1:length(gi)
                newlen = min(find(diff(good(run0:end,gi(i)))));
                if newlen > bestlen
                    bestdir = gi(i);
                    bestlen = newlen;
                end
            end
            if bestlen > 1
                % Now we can save this run of bestlen points on bestdir
                % Do we save the turnaround point? Why not? If it isn't
                run1 = run0 + bestlen - 1;
                newdirs(isbus(run0:run1)) = bestdir;
                n_runs = n_runs+1;
                runs(n_runs) = struct('busID', buses(bus), 'dir', bestdir, ...
                    'time', A.time(isbus(run0:run1)), ...
                    'dist', dist(run0:run1,bestdir), ...
                    'xdist', xdist(run0:run1,bestdir));
                run0 = run1;
            else
                run0 = run0+1;
                fprintf(1,'%s: Discarding single point\n', buses{bus});
            end
        end
    end
end
%% Graph runs by direction
for dir = 1:length(route.direction)
    figure;
    for run = find([runs.dir] == dir)
        plot(runs(run).time,runs(run).dist,'.-');
        hold on;
    end
    hold off;
    title([route.route ': ' route.direction(dir).title ' ' route.direction(dir).name ]);
end

%% Now routes are cleared up, but we may still need to make some
% breaks if a bus backtracks.
% showing all the buses for a particular direction
% instead of splitting on gaps, let's split whenever a bus
% appears to backtrack on the route.
% for dir = 1:length(route.direction)
%     P = route.direction(dir).stops(:,[2 3]);
%     dist = route.direction(dir).dist;
%     cdist = route.direction(dir).cdist;
%     figure;
%     for bus = 1:length(buses)
%         % find all the points where this bus is on this direction
%         isbus = find(strcmp(A.ID, buses(bus) ));
%         xi = isbus(dirs(isbus) == dir);
%         % then calculate distance along the direction
%         d = zeros(length(xi),1);
%         for j = 1:length(xi)
%             xdist = distance([A.lat(xi(j)) A.lon(xi(j))], P);
%             st = find(xdist == min(xdist));
%             % I'll assume only one.
%             st = st(1);
%             % Given P1, P2 on route and bus position B
%             % compare dist(P1,B)+dist(B,P2) with dist(P1,P2)
%             % (d1+d2-d0)/d0 and choose the one closer to zero
%             if st ~= 1 && st == length(xdist)
%                 st = st - 1;
%             elseif st ~= 1
%                 q0 = xdist(st)+xdist(st-1)-dist(st-1);
%                 q1 = xdist(st)+xdist(st+1)-dist(st);
%                 if q0 < q1
%                     st = st-1;
%                 end
%             end
%             d(j) = cdist(st) + xdist(st);
%         end
%         d = d*rfact;
%         % now split whenever the bus backtracks by more than .25 miles
%         backtracks = diff(d)< 0;
%         btstart = find(diff([ 0; backtracks ]) > 0);
%         btend = find(diff([backtracks; 0]) < 0)+1;
%         tot_bt = d(btstart)-d(btend);
%         bts = find(tot_bt >= .25);
%         n_bts = length(bts);
%         cumtime = NaN*zeros(length(xi)+n_bts);
%         cumdist = NaN*zeros(length(xi)+n_bts);
%         cumi = 1;
%         for i = 0:n_bts
%             if i == 0
%                 if n_bts == 0
%                     xxi = 1:length(xi);
%                 elseif bts(1) ~= 1 || btstart(1) ~= 1
%                     xxi = 1:btstart(bts(1))-1;
%                 else
%                     xxi = [];
%                 end
%             elseif i < n_bts
%                 xxi = btend(bts(i)):btstart(bts(i+1));
%             elseif btend(bts(i)) < length(xi)
%                 xxi = btend(bts(i)):length(xi);
%             else
%                 xxi = [];
%             end
%             if ~isempty(xxi)
%                 newcumi = cumi+length(xxi);
%                 cumtime(cumi:newcumi-1) = A.time(xi(xxi));
%                 cumdist(cumi:newcumi-1) = d(xxi);
%                 cumi = newcumi+1;
%             end
%             % h = plot(A.time(xi), d);
%         end
%         h = plot(cumtime, cumdist);
%         set(h, 'DisplayName', buses{bus} );
%         hold on;
%     end
%     hold off;
%     title([route.route ': ' route.direction(dir).title ' ' route.direction(dir).name ]);
% end
