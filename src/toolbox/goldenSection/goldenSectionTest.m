function goldenSectionTest(  )
%goldenSectionTest tests goldenSection to verify its correct implementation
xtol = 1e-4;
[x,fval,debug] = goldenSection (@(x) costFun_convex1D (x),-2,4.5,xtol,[],[costFun_convex1D(0) costFun_convex1D(4)]);
fprintf('Function Cost = %0.8e \n',fval);
fprintf('Optimum X     = %0.8e \n',x);
fprintf('Error   X     = %0.8e \n',abs(x-2));
fprintf('Bracket Distance = %0.8e \n',debug.bracketLength);
fprintf('Number Of Iterations  = %i \n',debug.numberOfIterations);
fprintf('Number Of Evaluations = %i \n',debug.numberOfEvaluations);
end

