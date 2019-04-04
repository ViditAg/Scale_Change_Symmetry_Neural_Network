function [ ava cl ] = avalanche(pp)
%LOCALAVALANCHEFINDER Compute cluster and avalanche data from a point process
% pp is a (x,y,t) point process timecourse
% This function returns a description of clusters (cl)
% and avalanches (av) detected in pp
% cl.C is a cell array of clusters (connected components) at each t
% cl.AV is a cell array of active voxel ids at each t
% cl.Lab is a 4D timecourse where (x,y,z,t) is a cluster label
% cl.N is a vector of the number of clusters at each t
% cl.A is a vector of the number of active voxels at each t
% cl.P is a vector of the order parameter (size of the largest cluster) at
% each t
% cl.Pmag is a vector of the order parameter (size of the largest cluster)
% at each t based on MAGNITUDE of the values within it
% cl.S is a vector of the frequency of cluster sizes 1..(x*y*z)
% cl.D is a cell array of the fractal dimensions of clusters at each t
%   ** D is EXPERIMENTAL **
% ava.N is the number of avalanches detected in total
% ava.O is a vector of time of onset of each avalanche
% ava.L is a vector of the duration of each avalanche
% ava.S is a vector of the size of each avalanche
% ava.St is a cell array of the size of each avalanche over time
% ava.M is a vector of the maximum size
% ava.A{i}{t} is a cell array of the voxel ids belonging to each avalanche i at
% each time point for that avalanche
%
%
%   by Gregory Scott (gregory.scott99@imperial.ac.uk)
%   based on Tagliazucchi et al, Frontiers in Physiology, 2012

%-------------------------------------------------------------------------
% Cluster analysis

cl.ImageSize = [ size(pp,1) size(pp,2) ];
cl.S = zeros(size(pp,1) * size(pp,2),1); % create a frequency table for cluster sizes in VOXELS

clusterim = zeros(size(pp(:,:,1))); % create a volume for box counting

for t=1:size(pp, 3) % iterate over pp timecourse
    im = squeeze(pp(:,:,t)); % pull out image at this time point
    ConComp = bwconncomp(im); % find connected components
    
    cl.C{t} = ConComp; % store the clusters (connected components)
    cl.N(t) = ConComp.NumObjects; % number of clusters
    cl.A(t) = nnz(im); % number of active sites
    cl.AV{t} = vertcat(ConComp.PixelIdxList{:}); % list of voxel ids
    %cl.P(t) = 0; % size of largest cluster by VOXELS (start at zero)
    %cl.Pmag(t) = 0; % size of largest cluster by AMPLITUDES
    %cl.D{t} = []; % fractal dimensions
    labelim = zeros(size(pp(:,:,1))); % create a volume for labelling
    
    % iterate over each cluster, recording the statistics
    for i=1:ConComp.NumObjects
        sz = length(ConComp.PixelIdxList{i}); % size of cluster in VOXELS
        
        % calculate size of cluster in MAGNITUDE
        clMag{t,i} = sum(im(ConComp.PixelIdxList{i}));
        cl.S(sz) = cl.S(sz) + 1; % update frequency counts for size
        %cl.P(t) = max(cl.P(t), sz); % update size of largest cluster by VOXELS
        %cl.Pmag(t) = max(cl.Pmag(t), clMag{t,i}); % update by MAGNITUDE
        labelim(ConComp.PixelIdxList{i}) = i; % light up labels
%         if(sz == 1) % box counting unnecessary if single voxel
%             cl.D{t} = [ cl.D{t} 2 ]; % assume dimension = 3 (?)
%         else % box counting and fractal dimension
%             % light up voxels for this cluster in a temporary image
%             clusterim(ConComp.PixelIdxList{i}) = 1;
%            % fd = boxcounterik(trimmask(clusterim)); % box count
%             clusterim(ConComp.PixelIdxList{i}) = 0; % turn off the cluster
%          %   cl.D{t} = [ cl.D{t} fd ]; % store the fractal dimension
%         end
    end
    cl.Lab(:,:,t) = labelim; % store cluster labelling for this time point
end

%-------------------------------------------------------------------------
% Avalanche analysis v2 - working backwards

%     Let Cti be the ith cluster at time t. We consider a cluster i0 starting
%     an avalanche at time t0 if for all j,
%
%     Ct0?1j n Ct0i0 = 0 (i.e., no clusters were present in that region at the previous timestep)
%
%     An id is assigned to this avalanche and the same id is assigned to all clusters
%     intersecting this cluster at the following time,this is all clusters i such that
%
%     Ct0i0 ? Ct0+1i != ?.
%
%     The same procedure is applied recursively to all clusters satisfy-
%     ing the former condition until no more intersections are found.
%
%     When this happens, all clusters labeled with this id constitute the
%     avalanche.

ava.N = 0; % avalanche counter (and identifier)
%tic
for t=2:size(pp,3)
   % toc
    for i = 1:cl.C{t}.NumObjects % iterate over clusters at this time point
        
        % is this cluster the start of a new avalanche? to be one,
        % all the voxels in the cluster
        % must have been off in the previous time point
        % AND ALL VOXELS ADJACENT TO THE CLUSTER (ELSE THEY WOULD HAVE BEEN
        % IN THE AVALANCHE AND SO NOT MARK THE START OF A NEW AVA)
        
        % Get the voxels in the cluster
        voxels = cl.C{t}.PixelIdxList{i};
        
        % dilate the cluster by one voxel
        voxels2 = dilator2(voxels, cl.ImageSize);
        % Test whether the cluster intersects with
        % no active voxels at the previous time point
        % TO DO: might be faster just to do something like
        % img=pp(:,:,t-1); img(voxels) == 0;
        if(isempty(fastintersect(voxels2, cl.AV{t-1})))
            % this must be a new avalanche to track
            ava.N = ava.N + 1; % increment avalanche counter
            id = ava.N;
            ava.O(id) = t; % record time of onset
            ava.S(id) = length(voxels); % record number of VOXELS
            %ava.Smag(id) = clMag{t,i}; % record size in MAGNITUDE
            %ava.St{id}(1) = length(voxels);
            %ava.L(id) = 1; % record duration (start at 1)
            % TO DO: uncomment this if ever required (slow!)
          %  ava.A{id, 1} = voxels; % record ids of voxels at avalanche onset
           % ava.Voxels{id} = voxels;
            
            % track this avalanche through time and see which clusters
            % intersect it and when it comes to an end
            for t2 = (t+1):size(pp,3)
                
                % see if there is any intersection of the
                % voxels of the avalanche from the previous time point
                % with any active voxels in the present
                
                % dilate the avalanche by one voxel before checking for
                % intersections
                voxels = dilator2(voxels, cl.ImageSize);
                
                % use the labelling lookup volume to find intersecting clusters
                labelim = squeeze(cl.Lab(:,:,t2));
                intersectingids = labelim(voxels);
                
                % are any labels non-zero (i.e. clusters)?
                if(any(intersectingids))
                    % the avalanche has survived this time point!
                    %ava.L(id) = t2 - t; % record new duration
                    
                    % remove zeros (shouldnt be any) and duplicate labels
                    intersectingids = unique(intersectingids(intersectingids > 0));
                    
                    % intersectingids is a list of unique clusters ids which
                    % have contiguity with the (dilated) avalanche, so we extend
                    % the avalanche to include all voxels in these clusters
                    % (i.e. the avalanche propagates)
                    newvoxels = vertcat(cl.C{t2}.PixelIdxList{intersectingids});
                    
                    % TO DO: uncomment this if ever required (slow!)
                   % ava.A{id, (t2-t) + 1} = newvoxels; % store the new voxels
                   % ava.Voxels{id} = [ ava.Voxels{id}; newvoxels ];
                    
                    % AVA SIZE USING NUMBER OF ACTIVE PIXELS IN AVA AT EVERY
                    % TIME POINT
                    ava.S(id) = ava.S(id) + length(newvoxels);
                    
                    % Calculate new size in magnitude
                    %ava.Smag(id) = ava.Smag(id) + sum([clMag{t2,intersectingids}]);
                    % store the progression in size over time
                    % TO DO: this is duplicating the role of ava.A
                   % ava.St{id}((t2-t) + 1) = length(newvoxels);
                    
                    voxels = newvoxels; % update the voxels for the next time point
                    
                else
                    %% NEW
                    % Calculate new ava S AVA SIZE USING UNIQUE PIXELS IN AVA AT ANY TIME POINT
                    %ava.S(id) = length(unique(ava.Voxels{id}));
                   %  ava.Voxels{id} = unique( ava.Voxels{id} );
                    break;
                end
            end
        end
    end
end
if ava.N == 0
   ava.O = []; 
   ava.S = []; 
  % ava.St = []; 
  % ava.L = []; 
   ava.Voxels = [];   
end
end