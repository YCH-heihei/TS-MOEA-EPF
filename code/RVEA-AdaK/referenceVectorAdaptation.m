function r = referenceVectorAdaptation(Population,zmin)
    [Fno,~]=NDSort(Population.objs,1);
    zmax=max(Population(Fno==1).objs,[],1); % only non-dominated solutions in the population.
    r=zmax-zmin;
end