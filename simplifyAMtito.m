function [nodes,HAM,DAM,OAM]=simplifyAM(nodes,HAM,DAM,OAM)
% Prunes out 2nd degree nodes and the intermediate edges from HAM, DAM
%
% DETAIL:
%           The script prunes out 2nd degree nodes and the similar intermediate
%           edges from HAM and DAM. Additionally, the distance of the
%           pruned edges are resummed in DAM. Only the edges with equal
%           weight in HAM are pruned.
% INPUT:
%           HAM(i,j) (Sparse) - Adjacency matrix containing the highway
%               class between nodes i and j
%           DAM(i,j) (Sparse) - Adjacency matrix containing the distance
%               information between nodes i and j
%           nodes (Double x 2) - Node longitude and latitude of the array
%               index for reference by HAM & DAM
% OUTPUT:
%           HAM(i,j) (Sparse) - Adjacency matrix containing the highway
%               class between nodes i and j
%           DAM(i,j) (Sparse) - Adjacency matrix containing the distance
%               information between nodes i and j
%           nodes (Double x 2) - Node longitude and latitude of the array
%               index for reference by HAM & DAM
% EXAMPLE:
%           [nodes,HAM,DAM,OAM] = removeAM2DegreeNodes(nodes,HAM,DAM,OAM)
% AUTHOR:
%           Bharat Kunwar
%           https://github.com/bkunwar/AMTools

if ~exist('OAM','var')
    OAM = false;
end

% Before state
beforeEdges=length(find(DAM));
beforeNodes=length(nodes);

% This is the first run
iteration = 0;

% Initiate a progress bar
step = 'Processing 2nd degree nodes...';
h = waitbar(0,step);

before = sum(sum(DAM));

% Continue if it is the first run or if the last run yielded some work!
while iteration == 0 || work > 0
    save(['state/' datestr(now)]);
    
    % Reset work counter to zero
    work = 0;

    % Count the number of iterations
    iteration = iteration+1;
    
    % Express the adjacency matrix in terms of 0s and 1s
    logicalHAM = logical(HAM);

    % Calculate the degree of the adjacency matrix
    iDegree = full(sum(logicalHAM));
    oDegree = full(sum(logicalHAM'));
    
    % Find all the nodes with two in degree and two out degree
    tito = find(and(iDegree==2,oDegree==2));
    
    num_tito = length(tito);
    
    throughNodes = [];
    totalThisTurn = 0;
    
    %% Iterate over tito
    for i = 1:num_tito
        % Find the adjacent nodes connected to the node i on the left and right
        this_tito = tito(i);
        [pred,~,pred_value]=find(HAM(:,this_tito));
        [~,succ,succ_value]=find(HAM(this_tito,:));
        
        % These are our neighbours
        neig = unique([pred' succ]);
        
        % Only if the neighbours are not pending removal,
        % it is okay to remove this node
        unlistedNeighbours = isempty(intersect(neig,throughNodes));
        
        % Make sure the neighbours are the same types of edges according to HAM
        roadType = unique([pred_value' succ_value]);
        identicalNeighbours = length(roadType)==1;
        
        % The node pending removal can only have a neighbour on either side
        twoNeighbours = length(neig)==2;
        
        % Make sure that there aren't already preexisting connections
        % between the predecessors/successors
        unconnectedNeighbours = ~(DAM(pred(1),pred(2))||DAM(pred(2),pred(1)));
        
        % Ensure that the weight of the edges are the same
        if twoNeighbours && identicalNeighbours && unconnectedNeighbours %&&unlistedNeighbours
            distance = sum(DAM(:,this_tito));
            
            HAM(pred(1),pred(2))=roadType;
            HAM(pred(2),pred(1))=roadType;
            
            DAM(pred(1),pred(2))=distance;
            DAM(pred(2),pred(1))=distance;
            
            % Since this node meets all the required condition, flag it up to really remove it 
            throughNodes = [throughNodes this_tito];
            
            % Proof of work
            work = work + distance;
        end
        
        % Update the progress bar
        if mod(i,100) == 0 || mod(i,num_tito) == 0
            completed = i/num_tito;            
            progress = [num2str(i) ' of ' num2str(num_tito) ' complete.'];
            waitbar(completed,h,[step progress]);            
        end
    end
    
    % Remove the nodes pending removal
    HAM(:,throughNodes)=[];
    HAM(throughNodes,:)=[];            
    DAM(:,throughNodes)=[];
    DAM(throughNodes,:)=[];
    if issparse(OAM)
        OAM(:,throughNodes)=[];
        OAM(throughNodes,:)=[];    
    end
    nodes(throughNodes,:)=[];

    after = sum(sum(DAM));
    
    % Show percentage of the nodes that are of 2 degree    
    disp(['Iteration ' num2str(iteration) ': ' num2str(num_tito) ' of ' num2str(length(nodes)) ' nodes are currently tito nodes.']);
    disp([num2str((before-after)/after*100) '% of original road length has disappeared.']);
    
    % Show the distribution of node degrees before and after the simplification    
    if iteration == 1 || work == 0  
        figure;
        hist(iDegree+oDegree);
        figure;
        wgPlot(HAM,nodes,'vertexMarker','none');
    end
end

% Close the progress bar
close(h);

% After state
afterEdges=length(find(DAM));
afterNodes=length(nodes);

% Output the results
disp(['Edges - Before: ' num2str(beforeEdges) ' | After: ' num2str(afterEdges) ' | Total pruned: ' num2str(beforeEdges-afterEdges)]);
disp(['Nodes - Before: ' num2str(beforeNodes) ' | After: ' num2str(afterNodes) ' | Total pruned: ' num2str(beforeNodes-afterNodes)]);
disp([ num2str(num_tito) ' nodes cannot be removed.']);