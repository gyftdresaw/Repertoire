function plotData(X)
%data = dlmread('Compare_V_genes.txt', '\t')
%X = data([2:end], 2)
%PLOTDATA Plots the data points X into a new figure 
%   PLOTDATA(x) plots the data points with + for the positive examples
%   and o for the negative examples. X is assumed to be a Mx1 matrix.
SPO = zeros(5,1);
SPY = zeros(5,1);
THYO = zeros(5,1);
THYY = zeros(5,1);
% Create New Figure
figure; hold on;

% ====================== YOUR CODE HERE ======================
% Instructions: Plot the positive and negative examples on a
%               2D plot, using the option 'k+' for the positive
%               examples and 'ko' for the negative examples.
%


SPO = X([1:5]);
SPY = X([6:10]);
THYO = X([11:15]);
THYY = X([16:20]);

Vals  = [mean(SPO), mean(SPY), mean(THYO), mean(THYY)];
Err = [std(SPO), std(SPY), std(THYO), std(THYY)];

errorbar(Vals, Err)


% =========================================================================



hold off;

end
