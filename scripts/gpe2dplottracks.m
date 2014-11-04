function gpe2dplottracks(tracks,adjacency_tracks, points)

n_tracks = numel(tracks);

%for i = 1:length(points);
%    points{i} = [points{i},repmat(i,length(points{i}),1)];
%end
all_points = vertcat(points{:});


for i_track = 1 : n_tracks
    track = adjacency_tracks{i_track};
    track_points = all_points(track, :);
    if(track_points(1,3)<0 && length(track_points(:,1))>10)
        plot(smooth(track_points(:,1)',20,'moving'),smooth(track_points(:,2)',20,'moving'),'LineWidth',2);
    end
    
   if(track_points(1,3)>0 && length(track_points(:,1))>10)
        plot(smooth(track_points(:,1)',20,'moving'),smooth(track_points(:,2)',20,'moving'),'r','LineWidth',2);
   end
    
        
    hold on
end
end