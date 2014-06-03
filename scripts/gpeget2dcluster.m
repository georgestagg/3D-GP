function [clusters] = gpeget2dcluster(xlocs,ylocs,pol)
clusters{1} = [];

points =  [xlocs;ylocs]';
polarity = pol';
D = pdist2(points,points);
rmdipoles = 1;
clstring = 1;
while(rmdipoles)
    %get nn for each vortex
    nnList = knnsearch(points,points,'k',2);
    nnList = nnList(:,2);
    for i=1:length(nnList)
      if(i == nnList(nnList(i)) && polarity(i) ~= polarity(nnList(i)))
          %disp('Opposite sign mutal NN found, removing both')
          points([i,nnList(i)],:) = [];
          polarity([i,nnList(i)]) = [];
          break;
      end
      if(i == length(nnList))
        %we made it to the end with no more dipoles removed!
        %disp('Finished removing Dipoles')
        rmdipoles=0;
      end
    end
end

%get nn for each vortex
nnList = knnsearch(points,points,'k',2);
nnList = nnList(:,2);
D = pdist2(points,points);
for i=1:length(nnList)
    for j=1:length(nnList)
        if(i==j)break;end;
         if(polarity(i) == polarity(j))
             distij = D(i,j);
             distsi = D(i,:);
             dists_nPoli = distsi(polarity == -polarity(i));
             distsj = D(j,:);
             dists_nPolj = distsj(polarity == -polarity(j));
             minDistTonPol = min([dists_nPoli,dists_nPolj]);
             
             if(distij < minDistTonPol)
                %disp('Clustering two vortices...')
                clusters=clusterup(clusters,i,j);
             end
         end
    end
end

clusters = clusters(~cellfun('isempty',clusters));  
%for cl=1:length(clusters)
%    clusters{cl} =  points(clusters{cl},:);
%    plot(clusters{cl}(:,1),clusters{cl}(:,2))
%end


    function [clusters]=clusterup(clusters,a,b)
        founda = 0;
        foundb = 0;
        for k=1:length(clusters)
            if (founda == 0 && sum(clusters{k} == a)>0)
                founda = k;  
            end
            if (foundb == 0 && sum(clusters{k} == b)>0)
                foundb = k;
            end
        end
        if(founda >0 && foundb == 0)
            %fprintf('Added vortex %d to cluster %d\n',b,founda);
            clusters{founda} = [clusters{founda},b];         
        end
        if(founda == 0 && foundb > 0)
            %fprintf('Added vortex %d to cluster %d\n',a,foundb);
            clusters{foundb} = [clusters{foundb},a];      
        end       
        
        if(founda > 0 && foundb > 0 && founda ~=foundb)
            %fprintf('Joining cluster %d and %d\n',founda,foundb,a,b);
            clusters{foundb} = [clusters{foundb},clusters{founda}];      
            clusters{founda} = [];
        end        
        
        if (founda == 0 && foundb == 0)
            %fprintf('Added vortex %d and %d to a new cluster, %d\n',a,b,length(clusters)+1);
            clusters{length(clusters)+1} = [a,b];
        end
    end
end

