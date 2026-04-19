function [A,Wadapt,P,B] = AdaK(A, Q, P, Wadapt, W, N, zmin, fstep, flast, generations, max_generations, flag, znadir, scale, update, T, B, utility)
    if(generations<flast*max_generations)
        if(length(znadir)==0)
            % 非被支配个体构成的极差点。
            % [Fno,~]=NDSort(A,1);
            [Fno,~]=NDSort(P.objs,1);
            zmax=max(P(Fno==1).objs,[],1);
        else
            zmax=znadir;
        end
        zMean=mean(A,1);
        zStd=std(A,1);
        zppf=zMean+6*zStd;

        [A,S,Wvalid]=nichingPairPotentialArchive(A, Q, W, N, zmin, zmax, zppf, flag);
        if(mod(generations,ceil(fstep*max_generations))==0)
           % 需要对参考向量进行更新。
           Wadapt = adaptReferenceSet(A, S, Wvalid, W, N, zmin, zmax, scale);
           if update
              % 需要对小生境进行更新。
              [P, B] = updateBestAndNeighbors(P, Wadapt, N, T, zmin, zmax, scale, utility);
           
           else
              P=[];
              B=[];
           end
       end
    end
end