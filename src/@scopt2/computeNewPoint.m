function Z = computeNewPoint(~,searchDirection,stepLength,Z_p)
% Computes new point along search direction
Z  = Z_p + searchDirection*stepLength;
end

