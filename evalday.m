function B = evalday( rt, day )
route = eval(['route' rt]);
% coslat = cosd(P(:,1));
% CP = [ coslat .* cosd(P(:,2))  coslat .* sind(P(:,2))  sind(P(:,1)) ];
% N = cross(CP(1:end-1,:),CP(2:end,:));
% N = diag(1./sqrt(sum(N.*N,2))) * N;
% D = cross(CP(1:end-1,:),N); % should be pre-normalized
% dist2 = abs(asind(dot(CP(2:end,:),D,2)));
% x = 1:size(N,1);
% plot(x,dist,'+',x,dist2,'*');
A = ReadBusLog(['route_' rt '/' rt '_' day '.log']);
A.time = time2d(A.time)/3600-5;
% Fill in direction information
dirs = zeros(size(A.time));
for dir = 1:length(route.direction)
    dirs(strcmp(route.direction(dir).tag,A.dirTag)) = dir;
    P = route.direction(dir).stops(:,[2 3]);
    route.direction(dir).dist = distance(P(1:end-1,:),P(2:end,:));
    route.direction(dir).cdist = [ 0; cumsum(route.direction(dir).dist) ];
end
% Unknown directions have index of zero
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
% Parse bus route into runs and assign directions
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
n_runs = 0;
runsidx = [];
B.dist = [];
B.xdist = [];
B.run = [];
B.newdirs = [];
newdirs = dirs*0;
for bus = 1:length(buses)
    isbus = find(strcmp(A.ID, buses(bus) ));
    [ dist, xdist ] = evaldir(A, route, isbus);
    fwd = diff(dist,1,1) > -.01;
    if isempty(fwd)
        fwd = zeros(1,length(dist));
        fwdfrom = fwd;
    else
        tgap = (diff(A.time(isbus),1,1)>.2);
        fwd(tgap,:) = 0;
        % if the direction is right and the time gap is small enough,
        fwdfrom = [fwd; fwd(end,:)];
        fwd = [ fwd(1,:); fwd ];
    end
    near = xdist < .02;
    good = fwd & near;
    goodfrom = fwdfrom & near;
    good = [ good; 0*good(end,:) ];
    run0 = 1;
    while run0 < length(isbus)
        % Find the next good point
        good(run0,:) = goodfrom(run0,:);
        run1 = find(any(good(run0:end,:)'), 1 );
        if isempty(run1)
            fprintf(1,'%s: Discarding %d points at end of run\n', ...
                buses{bus}, length(isbus)-run0+1);
            A.dist(isbus(run0:end)) = dist(run0:end,1);
            A.xdist(isbus(run0:end)) = xdist(run0:end,1);
            run0 = length(isbus);
        else
            if run1 > 1
                % fprintf(1, '%s: Discarding %d points\n', buses{bus}, run1-1);
                A.dist(isbus(run0:run1-1)) = dist(run0:run1-1,1);
                A.xdist(isbus(run0:run1-1)) = xdist(run0:run1-1,1);
                run0 = run0 + run1 - 1;
            end
            % Now select the longest run
            gi = find(good(run0,:));
            bestdir = 0;
            bestlen = 1;
            for i = 1:length(gi)
                newlen = find(diff(good(run0:end,gi(i))), 1 );
                if newlen > bestlen
                    bestdir = gi(i);
                    bestlen = newlen;
                end
            end
            if bestlen > 1
                % Now we can save this run of bestlen points on bestdir
                % Do we save the turnaround point? Why not?
                % Problem is when it gets assigned to a second run,
                % it gets a second dist, which is probably near zero
                % To do this right, I should collect the runs, and
                % then reconstruct the points, duplicating any
                % turnarounds.
                run1 = run0 + bestlen - 1;
                % newdirs(isbus(run0:run1)) = bestdir;
                n_runs = n_runs+1;
                newdirs(isbus(run0:run1)) = bestdir;
                runsidx = [ runsidx; isbus(run0:run1) ];
                B.dist = [B.dist; dist(run0:run1,bestdir)];
                B.xdist = [B.xdist; xdist(run0:run1,bestdir)];
                B.newdirs = [ B.newdirs; ones(run1-run0+1,1)*bestdir];
                B.run = [ B.run; ones(run1-run0+1,1)*n_runs];
                run0 = run1;
            else
                A.dist(isbus(run0)) = dist(run0,1);
                A.xdist(isbus(run0)) = xdist(run0,1);
                fprintf(1,'%s: Discarding single point\n', buses{bus});
                run0 = run0+1;
            end
        end
    end
end
% Now should collect all the points that were not assigned to
% a direction.
B.ID = A.ID(runsidx);
B.dirTag = A.dirTag(runsidx);
B.time = A.time(runsidx);
B.age = A.age(runsidx);
B.lat = A.lat(runsidx);
B.lon = A.lon(runsidx);
B.heading = A.heading(runsidx);
B.dirs = dirs(runsidx);
B.route = route;
