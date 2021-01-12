function hfig =  plotSolutionMatrices(obj)
% Plots nonempty entries in solution matrices
Gaug = obj.problem.transcription.Gaug;
haug = obj.problem.transcription.haug;
MEQ  = obj.problem.transcription.MEQ;
pEQ  = obj.problem.transcription.pEQ;
f    = obj.problem.transcription.f;
fnp  = obj.problem.transcription.fnp;
dims = obj.problem.transcription.dims;

MINEQ = Gaug(1:dims.l,:);
pINEQ = haug(1:dims.l,1);
Gsoc  = Gaug(dims.l+1:end,:);
hsoc  = haug(dims.l+1:end,1);


scopt2.plotNonZeroEntriesAB(MEQ,pEQ,11,obj.problem.phases(1).index)
title('Equality Constraints')
scoptSave2pdf(obj,'equalityMatrices')
scopt2.plotNonZeroEntriesAB(MINEQ,pINEQ,12,obj.problem.phases(1).index)
title('Inequality Constraints')
scoptSave2pdf(obj,'inequalityMatrices')
scopt2.plotNonZeroEntriesAB(Gsoc,hsoc,13,obj.problem.phases(1).index)
title('Second Order Cone Constraints')
scoptSave2pdf(obj,'socMatrices')
scopt2.plotNonZeroEntries(f,14,obj.problem.phases(1).index)
title('Objective Vector')
scoptSave2pdf(obj,'objectiveVector')
scopt2.plotNonZeroEntries(fnp,14,obj.problem.phases(1).index)
title('Augmented Objective Vector')
scoptSave2pdf(obj,'augmentedObjectiveVector')

drawnow

hfig = gcf;

end


