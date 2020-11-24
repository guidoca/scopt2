function [x,fval,debug] = goldenSection (fhandle,a,b,xtol,numberOfIterations,Mbracket)
% Golden section search routine. It is a one dimensional optimiser which
% does not use derivative information and which can find the solution in a
% pre-determined number of iterations dependent on the variable tolerance.
% It searches within the inverval specified by a and b, and is based on the 'magic' ratio 0.6180
%
% fhandle              IN: Handle object of Objective function to be minimised
% a                    IN: lower interval boundary
% b                    IN: upper interval boundary
% xtol                 IN: variable tolerance for the search
% numberOfIterations   IN: number of iterations (if not inserted, computed
% with xtol)
% Mbracket             IN: Not necessary. Initial objective value at a and
% b, only if its desired to include boundaries in search
% x                   OUT: Optimum variable
% fval                OUT: Optimum objective
% debug               OUT: Additional debugging information as number of
% iterations
%
% Successive Convexification OPtimal ConTrol 2.0
% Optimal control tool to solve single phase autonomous problems
% using successive convexification.
% Copyright (C) 2020  Guillermo J. Dominguez Calabuig
% _________________________________________________________________________
%% Constants
goldenRatio     = (sqrt(5)-1)/2;
%% Read Input and Specifiy Algorithm settings
bracket       = [a,b]   ;
bracketLength = abs(b-a);
if ~exist('numberOfIterations','var')
    numberOfIterations = ceil(log(xtol / bracketLength) / log(goldenRatio));
elseif isempty(numberOfIterations)
    numberOfIterations = ceil(log(xtol / bracketLength) / log(goldenRatio));
end

if ~exist('Mbracket','var')
    Mbracket = [fhandle(a) , fhandle(b)];
    numberOfEvaluations =  2;
else
    if strcmp(Mbracket,'omit-boundaries')
        Mbracket = [nan nan];
    end
    numberOfEvaluations =  0;
end
%% Loop through steps
% Generate solution on both fractions to choose initial direction
x2   = bracket(1) + (1-goldenRatio)*(bracket(2)-bracket(1));
M2   = fhandle(x2) ;
x3   = bracket(1) + goldenRatio*(bracket(2)-bracket(1));
M3   = fhandle(x3) ;
numberOfEvaluations = numberOfEvaluations + 2;
% Iterate to find best triplet until bracketing distance is lower than xtol
for ii = 1:numberOfIterations
    if M2>M3 % chose minimum from triplets M2 M3 M4 and compute next M3
        bracket(1) = x2;
        Mbracket(1)= M2; 
        triplet  = [bracket(1) x3 bracket(2)];
        tripletM = [Mbracket(1) M3 Mbracket(2)];
        if ii <numberOfIterations
            x2  = x3;
            M2  = M3;
            x3  = bracket(1) + goldenRatio*(bracket(2)-bracket(1));
            M3  = fhandle(x3) ;numberOfEvaluations=numberOfEvaluations+1;
        end
    else % chose minimum from triplets M1 M2 M3 and compute next M2
        bracket(2) = x3;
        Mbracket(2)= M3;
        triplet  = [bracket(1) x2 bracket(2)];
        tripletM = [Mbracket(1) M2 Mbracket(2)];
        if ii < numberOfIterations
            x3  = x2;
            M3  = M2;
            x2  = bracket(1) + (1-goldenRatio)*(bracket(2)-bracket(1));
            M2  = fhandle(x2) ;numberOfEvaluations=numberOfEvaluations+1;
        end
    end
end
[fval,indAux] =   min(tripletM);
x = triplet(indAux);

%% Additonal output
debug.numberOfIterations  = numberOfIterations;
debug.numberOfEvaluations = numberOfEvaluations;
debug.bracketLength       = diff(bracket);
end

