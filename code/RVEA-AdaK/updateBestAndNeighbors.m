function  [Pupdate, B]=updateBestAndNeighbors(P, Wadapt, N, T, zmin, zmax, scale, type)
% Updates best solution and neighbors for each subproblem"""
    if scale
        denom = zmax-zmin;
        denom(denom == 0) = 1e-12;
        V = Wadapt./denom;
        denom = sum(V, [],1);
        denom(denom == 0) = 1e-12;
        Vweight = V./denom;
        bDist=squreform(pdist(Vweight));
        [~,orderIndex]=sort(bDis1,1);
        B=orderIndex(:,1:T);
    else
        bDist=squreform(pdist(Wadapt));
        [~,orderIndex]=sort(bDis1,1);
        B=orderIndex(:,1:T);
    end
    Pupdate = [];

    % MOEAD
    for i =1:N
        normW   = sqrt(sum(Wadapt(i,:).^2,2));
        normP   = sqrt(sum((P.objs-zmin).^2,2));
        % normO   = sqrt(sum((Offspring.obj-Z).^2,2));
        CosineP = sum((P.objs-zmin).*Wadapt(i,:),2)./normW./normP;
        % CosineO = sum(repmat(Offspring.obj-Z,T,1).*W(P,:),2)./normW./normO;
        fitness   = normP.*CosineP + 5*normP.*sqrt(1-CosineP.^2);
        
        I=randperm(1:size(P,1));
        [~,minIndex]=min(fitness(I));
        best = I(minIndex);
        Pupdate =[Pupdate;P(best)];
   end 
end