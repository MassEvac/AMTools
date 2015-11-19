fr = [1 2 3 4 4 4 5 6 6 7 8 9];
to = [2 1 2 3 5 6 4 4 7 6 7 8];

n = max([fr to]);

nodes = 1:n;
a = sparse(fr,to,1,n,n)

%      0     1     0     0     0     0     0     0     0
%      1     0     0     0     0     0     0     0     0
%      0     1     0     0     0     0     0     0     0
%      0     0     1     0     1     1     0     0     0
%      0     0     0     1     0     0     0     0     0
%      0     0     0     1     0     0     1     0     0
%      0     0     0     0     0     1     0     0     0
%      0     0     0     0     0     0     1     0     0
%      0     0     0     0     0     0     0     1     0

%[n2,a2,b2,c2] = simplifyAM([1:8; 1:8]',a,a);

view(biograph(a))

in = full(sum(a))
out = full(sum(a'))

% Node indices with one in one out
oioo = find(and(in==1,out==1))
% Node indices with two in two out
tito = find(and(in==2,out==2))



succ = @(DAM,n) find(adj(n,:));
pred = @(DAM,n) find(adj(:,n));
neig = @(DAM,n) unique([succ(adj,n) pred(adj,n)']);

removal = []

for i = oioo
    if succ(i) ~= pred(i)
        removal = [removal i];
    end
end



