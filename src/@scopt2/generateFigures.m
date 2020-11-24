function generateFigures(obj)
% Routine to generate output figures if sopecified by debugging
% and solution flags
obj.plotRealObjective();
obj.plotSolutionMatrices();
obj.plotAdaptiveTrustRegionMonitor();
obj.plotConvergenceHistory();

if obj.debugging.saveTrustRegion && obj.algorithm.sequential.activate
    obj.plotTrustRegionHistory();
    obj.plotTrustRegionErrorHistory();
end
if obj.debugging.saveVirtualControls && obj.algorithm.virtualControl.phases(1).include
    obj.plotVirtualControlHistory();
    obj.plotVirtualControlErrorHistory();
end
if obj.debugging.saveVirtualBuffers && obj.algorithm.virtualBuffer.phases(1).include
    obj.plotVirtualBuffersHistory();
    obj.plotVirtualBuffersErrorHistory();
end
if obj.debugging.generateFigureDefects
    obj.plotVerificationDefectConstraints(); % That defect constraints are zero (compute constraint satisfaction with virtual controls)
end
%             obj.plotVerificationPathConstraints(); % That buffer zones satisfy path constraints
obj.plotTime();
end

