%function [HAM,DAM,nodes]=removeAM2DegreeNodes(HAM,DAM,nodes)
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
%           [HAM,DAM,nodes] = removeAM2DegreeNodes(HAM,DAM,nodes)

% Before state
beforeEdges=length(find(DAM));
beforeNodes=length(nodes);

% This is the first run
iteration = 0;
figureNo = 1;

% Initiate a progress bar
step = 'Processing 2nd degree nodes...';
h = waitbar(0,step);

% Continue if it is the first run or if the last run yielded some work!
while iteration == 0 || work > 0
    % Reset work counter to zero
    work = 0;
    iteration = iteration+1;    
    
    % Express the adjacency matrix in terms of 0s and 1s
    logicalHAM = logical(HAM);

    % Calculate the degree of the adjacency matrix
    HAMsq=logicalHAM^2;
    nodeDegrees=full(diag(HAMsq));
    
    % Find all the nodes with 2 degrees which we want to remove
    nodesToRemove=find(nodeDegrees==2);

    % Show percentage of the nodes that are of 2 degree    
    disp(['Iteration ' num2str(iteration) ': ' num2str(length(nodesToRemove)) ' of ' num2str(length(nodes)) ' nodes are currently 2 degree nodes.']);

    numNodesToRemove = length(nodesToRemove);
    reallyRemove = [];
    %%
    for i = 1:numNodesToRemove
        % Find the adjacent nodes connected to the node i on the left and right
        thisNodeToRemove = nodesToRemove(i);
        [LHS,~,Lvalue]=find(HAM(:,thisNodeToRemove));
        [~,RHS,Rvalue]=find(HAM(thisNodeToRemove,:));

        % Ensure that the weight of the edges are the same
        if isequal(sort(Lvalue),sort(Rvalue'))&&length(LHS)==2
            thisNode = LHS(1);
            thatNode = RHS(RHS~=LHS(1));
            value = Lvalue(1);
            HAM(thisNode,thatNode)=value;
            HAM(thatNode,thisNode)=value;
            dist = full(sum(DAM(:,thisNodeToRemove)));
            DAM(thisNode,thatNode)=dist;
            DAM(thatNode,thisNode)=dist;
            
            % Since this node meets all the required condition, flag it up to really remove it 
            reallyRemove = [reallyRemove thisNodeToRemove];
            
            % Proof of work
            work = work + 1;            
        else
            nodesToRemove(i)=0;
        end

        if mod(i,100) == 0 || mod(i,numNodesToRemove) == 0
            completed = i/numNodesToRemove;            
            progress = [num2str(i) ' of ' num2str(numNodesToRemove) ' complete.'];
            waitbar(completed,h,[step progress]);            
        end
    end
    
    % Remove the nodes pending removal
    HAM(:,reallyRemove)=[];
    HAM(reallyRemove,:)=[];            
    DAM(:,reallyRemove)=[];
    DAM(reallyRemove,:)=[];
    nodes(reallyRemove,:)=[];

    % Show the distribution of node degrees before and after the simplification
    if iteration == 0 || work == 0    
        subplot(2,2,figureNo);
        figureNo=figureNo+1;        
        hist(nodeDegrees);
        subplot(2,2,figureNo);
        figureNo=figureNo+1;           
        wgPlot(HAM,nodes,'vertexMarker','none');
    end
    
    % Count the number of iterations
end

% Close the progress bar
close(h);



% After state
afterEdges=length(find(DAM));
afterNodes=length(nodes);

% Output the results
disp(['Edges - Before: ' num2str(beforeEdges) ' | After: ' num2str(afterEdges) ' | Total pruned: ' num2str(beforeEdges-afterEdges)]);
disp(['Nodes - Before: ' num2str(beforeNodes) ' | After: ' num2str(afterNodes) ' | Total pruned: ' num2str(beforeNodes-afterNodes)]);
disp([ num2str(numNodesToRemove) ' nodes cannot be removed.']);