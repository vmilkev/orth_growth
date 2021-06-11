% clf;
% figure(1);
% hold on;

R = 5;
max_time = 200;

cluster = {};

% instance of a first cell
%                     origin  radius
cluster{1,1} = mcell( [0 0 0], R );

deltaT = 1.0; % time increment

for t = 1:max_time
    
    sizeCluster = size( cluster,1 );
    
    % growth loop
    for i = 1:sizeCluster
        cluster{i,1}.growth( cluster,deltaT );
    end
    
    % if a cell is mature than it initiates a new cell
    % and starts elongating
    l = 1;
    for i = 1:sizeCluster
        if ( cluster{i,1}.isElongating && ~cluster{i,1}.isLinked )
            % add new cell to the cluster
            cluster{sizeCluster+l,1} = cluster{i,1}.addcell;
            l = l + 1;
        end
    end
    
    % visualize cluster growth: directions(confinements, growth), cells
    % and some data
    allcells = size( cluster,1 );
    if (allcells >= 1)
        
        clf;
        figure(1);
        hold on;
        axis equal
        for ii = 1:sizeCluster
            for l = 1:size( cluster{ii,1}.confinements,1 )
                cluster{ii,1}.vplot_conf(l);
            end
            cluster{ii,1}.vplot_dir(1.0);
            if ( cluster{ii,1}.conf_num < 3 )
                cluster{ii,1}.vplot_dir(-1.0);
            end
            %cluster{ii,1}.splot;
        end
        hold off;
        
        out = zeros(allcells, 8);
        for jj = 1:allcells
            out(jj,:) = [allcells cluster{jj,1}.name cluster{jj,1}.age cluster{jj,1}.conf_num cluster{jj,1}.isElongating cluster{jj,1}.isLinked cluster{jj,1}.currentElongation  cluster{jj,1}.isVirtual];
        end
        display(['         ' 'cells' '        ' 'id' '          ' 'age' '          ' 'conf' '        ' 'elong' '       ' 'linked' '      ' 'length' '      ' 'vitual'])
        format shortG
        out
    end
    
end

