chdir ('~/Dropbox/OSM')
load('cache/highway/osm_gb/Carlisle/DAM.mat')
load('cache/highway/osm_gb/Carlisle/HAM.mat')
load('cache/highway/osm_gb/Carlisle/OAM.mat')
load('cache/highway/osm_gb/Carlisle/nodes.mat')

before = full(sum(sum(DAM)));

[nodes2,HAM2,DAM2,OAM2]=simplifyAM(nodes,HAM,DAM,OAM);

after = full(sum(sum(DAM2)));

change = (before - after) / after * 100;