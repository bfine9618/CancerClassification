close all;
clc;
load('geneNetwork_rawPCNCI.mat');
load('signal_mutation.mat');
load('histology_subtype.mat');
gN = +gN;
%%%%%I renamed the gene network to gN%%%%%
%%%%%I renamed the signal mutation to sigMut%%%%%
%%%%%I renamed the histology_subtype to hist%%%%%

%% 2.1
%%%%%I renamed the gene network to gN%%%%%

loops = sum(diag(gN) ~= 0);
directed = ~issymmetric(gN);
wieghted = (max(max(gN)) > 1);

spy(gN);
title('Visualized Gene Network');

%%% With self loops
D = diag(sum(gN~=0,2));
lp = D - gN;

%%% Remove diag 
gNd = gN - diag(diag(gN));
g = graph(gNd);
l = laplacian(g);

%% 2.2
S = lp;

[V, D]  = eig(S);
d=diag(D);
[~,order]=sort(abs(d),'ascend');
V=V(:,order); d=d(order); D=diag(d);

if max(max(abs(S-V*D*V')))<1e-9
    S=V*D*V';
else
    error('The eigendecomposition is not good.\n')
end

variation=diag(V'*S*V);


figure;
plot(diag(D), variation); 
title('Eigenvalues vs Variation');
xlabel('$\lambda_{k}$','Interpreter','LaTeX');
ylabel('$TV$','Interpreter','LaTeX');

%%
D = diag(sum(gN));
lp = D - gN;

[V, ~] = eigs(lp, length(lp));
H = V';

gft = zeros(240,2458);

for i=1:240
    gft(i,:) = sigMut(i,:)*H;
end

type1 = gft(hist==1, :);
mu1 = mean(type1, 1);

type2 = gft(hist==2, :);
mu2 = mean(type2, 1);

n = sum(abs(gft), 1);

res = abs(mu1-mu2);
res = res./n;

figure;
plot(1:length(sigMut), abs(res));
title('Residual for each gene');
xlabel('Gene k');
ylabel('Residual');

figure;
boxplot(abs(res));
title('Residual for each gene');

%%

[k, acc] = kNN(sigMut, hist, 100);

plot(k, acc);
title('Nearest Neighbor Accuracy');
xlabel('k nearest comparisons');
ylabel('Accuracy %');

%%
[~, k] = max(gft, [], 2);

H1 = zeros(240, 2458);

for i=1:240
    H1(i,k(i)) = 1;
end

fgft = H1.*gft;
fsig = fgft * V;

[k, acc] = kNN(fsig, hist, 100);

plot(k, acc);
title('Filtered Nearest Neighbor Accuracy');
xlabel('k nearest comparisons');
ylabel('Accuracy %');

%%
p = [20 40 50 60 80 99];
accP = [];

[~, srt] = sort(res, 'descend');

figure;

for i=1:length(p)
    per = round((p(i)/100) * length(res));
    fltr = zeros(1, 2458);
    fltr(res*1000000 >= srt(per)) = res(res*1000000 >= srt(per)); 
    
    fgft = fltr.*gft;
    fsig = fgft * V;
    
    [k, acc] = kNN(fsig, hist, 10);
    
    accP = [accP acc'];
    
    plot(k, acc);
    hold on;
end

legend('20th Percentile', '40th Percentile', '50th Percentile','60th Percentile', '80th Percentile', '99th Percentile'); 
xlabel('K nearest neighbor Comparison');
ylabel('accuracy %');
title('Percentile Filtered kNN');


