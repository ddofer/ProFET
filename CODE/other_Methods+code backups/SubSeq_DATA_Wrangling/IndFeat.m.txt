% IndFeat.m - routine for selecting independent features for 
% modeling problems
%
% code by Will Dwinnell (predictr@bellatlantic.net)
%
% Sig = IndFeat(X,Y)
%
% Calculate significance level, 'Sig' of real variables (in columns) from 
% matrix 'X', based on their ability to distinguish 2 categories in column
% vector 'Y' (typically, variables are kept when significance >= 2.0).
% For a continuous output variable, it is suggested that values be split
% into lower and upper 50%.
%
% Uses Weiss/Indurkhya 'independent features' significance testing
% method- keep in mind that this is not intended to be the final 
% variable selection, just a means of eliminating obviously uninteresting 
% features. See Weiss' and Indurkhya's book, "Predictive Data Mining".
%
% Example:
%
% X = randn(5e3,25);
% Y = double(sum(X(:,1:3)')' > 1.1);
% IndFeat(X,Y)
% figure
% stem(IndFeat(X,Y))
% grid on
% zoom on
%
% Last modified: Nov-14-2006

function Sig = IndFeat(X,Y)

% find the two (more to come) class "names"
UniqueClass = unique(Y);

% find indexes for both classes
ClassA = (Y == UniqueClass(1));
ClassB = (Y == UniqueClass(2));
nA     = sum(double(ClassA));
nB     = sum(double(ClassB));

% calculate significances
Sig = ...
    abs(mean(X(ClassA,:)) - mean(X(ClassB,:)))  ./  ...
    sqrt(var(X(ClassA,:)) ./ nA + var(X(ClassB,:)) ./ nB);


% EOF


