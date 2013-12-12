% This script is a test to see if the results yielded by simplifyAM.m are valid.
%
% DETAIL:
%           It works by calculating the graphshortestpath between any two
%           randomly selected nodes in the original and simplified adjacency 
%           matrices and compare the difference in results within a defined
%           tolerance. To pass, there should be no difference in results.
% AUTHOR:
%           Bharat Kunwar
%           https://github.com/bkunwar/AMTools
% NOTE: 
%           Test runs seem to confirm that it is a valid simplification
%
% Assumes that the following variables have been defined:
%     DAMsimple
%     DAM
%     nodesSimple
%     nodes

% Load the AM for the relevant city
%     place = 'Bristol';
%     folderPath = '~/Dropbox/Coding/MATLAB/OSM/cache/_highway/';
%     load([folderPath place '/DAM.mat']);
%     load([folderPath place '/HAM.mat']);
%     load([folderPath place '/nodes.mat']);

% Retrieve the simplified matrix
%     [HAMsimple,DAMsimple,nodesSimple]=simplifyAM(HAM,DAM,nodes);

valid = 0;
tests = 100;

% Create an array of arbitrary nodes
nodeRange = [1 length(nodesSimple)];
originNodeSimple = randi(nodeRange, 1, tests);
destinNodeSimple = randi(nodeRange, 1, tests);

for i = 1:tests
    % Show the test number
    i

    % What are the coordinates in the nodesSimple 
    origin=nodesSimple(originNodeSimple(i),:);
    destin=nodesSimple(originNodeSimple(i),:);

    % What is the node number for the same coordinates in the original nodes
    originNode = find(ismember(nodes,origin,'rows'));
    destinNode = find(ismember(nodes,destin,'rows'));

    % Calculate the distance for original and simplified matrices
    actualDistance = graphshortestpath(DAM,originNode,destinNode);
    simpleDistance = graphshortestpath(DAMsimple,originNodeSimple(i),destinNodeSimple(i));

    % The output should yield the same results
    if (actualDistance - simpleDistance) < 0.001 % This is our tolerance
        valid = valid + 1;
    end
end

disp(['The simplification is valid ' num2str(valid) ' out of ' num2str(tests) ' times.']);