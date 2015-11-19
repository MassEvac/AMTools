function [nodes,HAM,DAM,OAM]=simplifyAM(nodes,HAM,DAM,OAM,place)
% Prunes out flow through nodes from HAM, DAM, OAM
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
%           [nodes,HAM,DAM,OAM] = simplifyAM(nodes,HAM,DAM,OAM)
% AUTHOR:
%           Bharat Kunwar
%           https://github.com/bkunwar/AMTools

if ~exist('OAM','var')
    OAM = false;
end

if ~exist('place','var')
    place = 'simplifyAM';
end

% Before state
beforeEdges=length(find(DAM));
beforeNodes=length(nodes);

% This is the first run
iteration = 0;

before = sum(sum(DAM));

% Continue if it is the first run or if the last run yielded some work!
while iteration == 0 || work > 0
    % Reset work counter to zero
    work = 0;

    % Count the number of iterations
    iteration = iteration+1;
    
    %% Express the adjacency matrix in terms of 0s and 1s
    logicalHAM = logical(HAM);

    % Calculate the degree of the adjacency matrix
    iDegree = full(sum(logicalHAM));
    oDegree = full(sum(logicalHAM'));

    if iteration == 1
        figure;
        hist(iDegree+oDegree);
        figure;
        wgPlot(HAM,nodes,'vertexMarker','none');
    end    
    
    throughNodes = [];
    
    %% Find all the nodes with two in degree and two out degree
    step = 'Processing 2in2out degree nodes...';
    tito = find(and(iDegree==2,oDegree==2));
    
    num_tito = length(tito);
    
    save(['state/' datestr(now)]);
    
    % Iterate over tito
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
        
        % Ensure that the weight of the edges are the same
        if twoNeighbours && identicalNeighbours &&unlistedNeighbours
            % Make sure that there aren't already preexisting connections
            % between the predecessors/successors
            unconnectedNeighbours = ~(DAM(neig(1),neig(2))||DAM(neig(2),neig(1)));
            if unconnectedNeighbours
                distance = sum(DAM(:,this_tito));
                antidist = sum(DAM(this_tito,:));
                
                disagreement = abs(distance - antidist);
                
                equalNeighbours = disagreement < 1e-09 ;
                % Where there are two nodes connected to the node in question
                % that are connected to the other two in both
                % directions, sometimes, a one way road may have collapsed
                % and the sum of edge lengths may be different. If that is
                % the case, we do not remove these nodes.
                if equalNeighbours
                    HAM(neig(1),neig(2))=roadType;
                    HAM(neig(2),neig(1))=roadType;
                    DAM(neig(1),neig(2))=distance;
                    DAM(neig(2),neig(1))=distance;
                    
                    % Since this node meets all the required condition, flag it up to really remove it 
                    throughNodes = [throughNodes this_tito];

                    % Proof of work
                    work = work + distance*2;
                else
                    disp([place ': Distance disagreement > 1e-09 (' num2str(disagreement) ') at node ' num2str(this_tito) '...']);
                end
            end
        end
        
        % Update the progress bar
        if mod(i,100) == 0 || mod(i,num_tito) == 0
            completed = i/num_tito;            
            progress = [num2str(i) ' of ' num2str(num_tito) ' complete.'];
            disp([place ': ' step progress]);
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
    
    %% Express the adjacency matrix in terms of 0s and 1s
    logicalHAM = logical(HAM);

    % Calculate the degree of the adjacency matrix
    iDegree = full(sum(logicalHAM));
    oDegree = full(sum(logicalHAM'));
    
    throughNodes = [];

    % [preds,~,preds_value] = find(HAM(:,oioo));
    % [~,succs,succs_value] = find(HAM(oioo,:));    
    
    %% Find all the nodes with one in degree and one out degree
    step = 'Processing 1in1out degree nodes...';
    oioo = find(and(iDegree==1,oDegree==1));
        
    num_oioo = length(oioo);

    save(['state/' datestr(now)]);    
    
    % Iterate over oioo
    for i = 1:num_oioo
        % Find the adjacent nodes connected to the node i on the left and right
        this_oioo = oioo(i);
        [pred,~,pred_value]=find(HAM(:,this_oioo));
        [~,succ,succ_value]=find(HAM(this_oioo,:));
        
        % These are our neighbours
        neig = unique([pred' succ]);
        
        % Only if the neighbours are not pending removal,
        % it is okay to remove this node
        unlistedNeighbours = isempty(intersect(neig,throughNodes));
        
        % Make sure the neighbours are the same types of edges according to HAM
        roadType = unique([pred_value' succ_value]);
        identicalNeighbours = length(roadType)==1;
        
        % The node pending removal must be connected to two other nodes
        twoNeighbours = length(neig)==2;
                
        % Ensure that the weight of the edges are the same
        if twoNeighbours && identicalNeighbours && unlistedNeighbours
            % Make sure that there aren't already preexisting connections
            % between the predecessors/successors
            unconnectedNeighbours = ~DAM(pred,succ);
            if unconnectedNeighbours
                distance = DAM(pred,this_oioo)+DAM(this_oioo,succ);
                
                HAM(pred,succ)=roadType;
                DAM(pred,succ)=distance;
                
                DAM(pred,this_oioo)=0;
                DAM(this_oioo,succ)=0;
                
                % Since this node meets all the required condition, flag it up to really remove it 
                throughNodes = [throughNodes this_oioo];

                % Proof of work
                work = work + distance;
            end
        end
        
        % Update the progress bar
        if mod(i,100) == 0 || mod(i,num_oioo) == 0
            completed = i/num_oioo;            
            progress = [num2str(i) ' of ' num2str(num_oioo) ' complete.'];
            disp([place ': ' step progress]);
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

    disp([place ': ' num2str((before-after)/after*100) '% of original road length has disappeared.']);    
    
    % Show percentage of the nodes that are of 2 degree    
    disp([place ': Iteration ' num2str(iteration) ': ' num2str(num_tito) ' 2in2out / ' num2str(num_oioo) ' 1in1out nodes left of ' num2str(length(nodes)) ' nodes.']);
    
    if work == 0  
        figure;
        hist(iDegree+oDegree);
        figure;
        wgPlot(HAM,nodes,'vertexMarker','none');
    end
end



% After state
afterEdges=length(find(DAM));
afterNodes=length(nodes);

% Output the results
disp([place ': Edges - Before=' num2str(beforeEdges) ' | After=' num2str(afterEdges) ' | Pruned=' num2str(beforeEdges-afterEdges)]);
disp([place ': Nodes - Before=' num2str(beforeNodes) ' | After=' num2str(afterNodes) ' | Pruned=' num2str(beforeNodes-afterNodes)]);
disp([place ': ' num2str(num_tito) ' 2in2out / ' num2str(num_oioo) ' 1in1out nodes are not through nodes.']);
