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
dirs = zeros(size(A.time));
for dir = 1:length(route.direction)
    dirs(strcmp(route.direction(dir).tag,A.dirTag)) = dir;
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

%% Now directions are cleaned up, so let's generate a crude plot
% showing all the buses for a particular direction
for dir = 1:length(route.direction)
    P = route.direction(dir).stops(:,[2 3]);
    dist = distance(P(1:end-1,:),P(2:end,:));
    cdist = [ 0; cumsum(dist) ];
    figure;
    for bus = 1:length(buses)
        isbus = find(strcmp(A.ID, buses(bus) ));
        isdir = dirs(isbus) == dir;
        dirstart = find(diff([ 0; isdir ]) > 0);
        dirend = find(diff([isdir; 0]) < 0);
        cumtime = NaN*zeros(length(isdir)+length(dirstart)-1);
        cumdist = NaN*zeros(length(isdir)+length(dirstart)-1);
        cumi = 1;
        for i = 1:length(dirstart)
            xi = isbus(dirstart(i):dirend(i));
            d = zeros(length(xi),1);
            for j = 1:length(xi)
                xdist = distance([A.lat(xi(j)) A.lon(xi(j))], P);
                st = find(xdist == min(xdist));
                % I'll assume only one.
                st = st(1);
                % Given P1, P2 on route and bus position B
                % compare dist(P1,B)+dist(B,P2) with dist(P1,P2)
                % (d1+d2-d0)/d0 and choose the one closer to zero
                if st ~= 1 && st == length(xdist)
                    st = st - 1;
                elseif st ~= 1
                    q0 = xdist(st)+xdist(st-1)-dist(st-1);
                    q1 = xdist(st)+xdist(st+1)-dist(st);
                    if q0 < q1
                        st = st-1;
                    end
                end
                d(j) = cdist(st) + xdist(st);
            end
            % rdist = distance(A.lat(xi(1:end-1)),A.lon(xi(1:end-1)), ...
            %     A.lat(xi(2:end)),A.lon(xi(2:end)));
            % rdist = [ 0; cumsum(rdist) ];
            % h = plot(A.time(xi), rdist);
            newcumi = cumi+length(xi);
            cumtime(cumi:newcumi-1) = A.time(xi);
            cumdist(cumi:newcumi-1) = d;
            cumi = newcumi+1;
            % h = plot(A.time(xi), d);
        end
        h = plot(cumtime, cumdist*rfact);
        set(h, 'DisplayName', buses{bus} );
        hold on;
    end
    hold off;
    title([route.route ': ' route.direction(dir).title ' ' route.direction(dir).name ]);
end
