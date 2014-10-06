function [drifter,dist] = gpeget2dvort_countdrift(total,tracks,adjacency_tracks, points)

n_tracks = numel(tracks);

for i = 1:length(points);
    points{i} = [points{i},repmat(i,length(points{i}),1)];
end
all_points = vertcat(points{:});


%disconsti = 9.8;
disconsti = 9.7;
disconsto = 9.85;

drifter = zeros(length(points),1);
for i_track = 1 : n_tracks
    track = adjacency_tracks{i_track};
    track_points = all_points(track, :);
 
    dists = sqrt(track_points(1,1).^2+track_points(1,2).^2);
    diste = sqrt(track_points(end,1).^2+track_points(end,2).^2);
    dist{i_track} = [dists,diste];
    
    if(dists <  disconsti && diste >  disconsti)
        if(track_points(end,4) > 85 && abs(dists-diste) > 3 )
        for i = track_points(end,4):length(points);
            drifter(i) = drifter(i) + 1;
        end
        end
    end
    if(dists >  disconsto && diste <  disconsto)
        if(track_points(end,4) > 85 && abs(dists-diste) > 3 )
        for i = track_points(end,4):length(points);
            drifter(i) = drifter(i) - 1;
        end
        end
    end
    
end

%for i = 1:length(points);
%    if(max(total)<total(i)+drifter(i))
%        drifter(i) = max(total)-total(i);
%    end
%end

repmax = repmat(total(85),1801,1);
for i = 1:85
    drifter(i) = 0;
    repmax(i) = 0;
end
drifter(length(points))=drifter(length(points)-1);

%clf
%bar(((0:1800)*0.0005*1000)/(15*2*pi),repmax)
%hold on
%bar(((0:1800)*0.0005*1000)/(15*2*pi),tsmovavg(total(85)-drifter','s',15))
%bar(((0:1800)*0.0005*1000)/(15*2*pi),tsmovavg(total,'s',15))
%plot(((0:1800)*0.0005*1000)/(15*2*pi),repmax)
%plot(((0:1800)*0.0005*1000)/(15*2*pi),tsmovavg(total(85)-drifter','s',15))
%plot(((0:1800)*0.0005*1000)/(15*2*pi),tsmovavg(total,'s',15))
%axis([0.45, 8, 0, total(85)])

clf
drifter2 = drifter(1510)-drifter';
stuff = total(85)-drifter' - total;
stuff = max(stuff(1000:1510))-stuff;
%plot(((85:1510)*0.0005*1000)/(15*2*pi),tsmovavg(drifter2','s',15));
%hold on
%plot(((85:1510)*0.0005*1000)/(15*2*pi),tsmovavg(stuff(85:1600),'s',15));

stuff=tsmovavg(stuff(70:1510),'s',15);
drifter2=tsmovavg(drifter2(70:1510),'s',15);
curvefit(((84:1510)*0.0005*1000)/(15*2*pi),stuff(15:end));
%curvefit(((84:1510)*0.0005*1000)/(15*2*pi),drifter2(15:end));

fclose('all');
end