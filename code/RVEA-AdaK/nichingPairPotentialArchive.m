function [A, S, Wvalid] = nichingPairPotentialArchive(A, Q, W, N, zmin, zmax, zppf, flag)
    % Updates external archive based on niching and pair-potential selection
    A=[A;Q.objs];
    % unique = np.sort(np.unique(np.around(A, 6), return_index=True, axis=0)[1])
    A=round(A,6);
    [A,~,~]=unique(A, 'rows', 'stable');
    [Fno,~]=NDSort(A,1);
    A = A(Fno==1,:);
    S = [];
    Wvalid = [];
    if(size(A,1) >= N)
        denom = zmax-zmin;
        denom(denom == 0) = 1e-12;
        % 标准化
        Aprime = (A-zmin)./denom;
        
        % 基于与参考向量的垂直距离来将个体关联到参考向量上。
        normOfAprime=sum(Aprime.^2,2).^(1/2);
        
        directionOfAprime=Aprime./(sum(Aprime,2));
        directionOfW=W./(sum(W,2));
        % normOfVectors=sum(W.^2,2).^(1/2);

        sinDis=(1-(1-pdist2(directionOfAprime,directionOfW,'cosine')).^2).^(1/2);
        
        perpendicularDistance=normOfAprime.*sinDis;
        [d,pi]=min(perpendicularDistance,[],2);
        % pi, d = associate(Aprime, W)
        rho=zeros(1,size(W,1));
        for i=1:size(Aprime,1)
           rho(pi(i))=rho(pi(i))+1;
        end
        % rho = nicheCount(pi, W)
        [A, S] = nichingPairPotentialSelection(A,  N, rho, pi, d, zmin, zppf, flag);
        Wvalid = W(rho > 0,:);
    end
end