function [A,S] = nichingPairPotentialSelection(A, N, rho, pi, d, zmin, zppf, flag)
% Selects N solutions from A using niching and pair-potential selection
    Choose = false(1,size(A,1));
    J = find(rho>0);
    for j =1:length(J)
        I = find(pi == J(j));
        I=I(randperm(length(I)));
        [~,minIndex]=min(d(I));
        Choose(I(minIndex)) = true;
    end

    S = A(Choose,:);
    selected = size(S,1);
    if (selected == N)
        A = S;
    else
        A=A(~Choose,:);
        A=A(all(A<zppf,2),:);

        A=[S;A];
        
        zmax_arch = max(A, [],1);
        denom = zmax_arch-zmin;
        denom(denom == 0) = 1e-12;
        Aprime = (A-zmin)./denom;
        Diss = dissimilarityMatrix(Aprime, flag);
        Memo = sum(Diss,1);
        while (size(A,1) > N)
            C = abs(Memo);
            % I = np.arange(selected, len(C))
            % np.random.shuffle(I)
            % worst = I[np.argmax(C[I])]
            I=selected+1:size(A,1);
            I=I(randperm(length(I)));
            [~,worst]=max(C(I));
            % try
                Diss(I(worst),:)=[];
            % catch
            %     temp=0;
            % end
            Diss(:,I(worst))=[];
            % Diss = np.delete(Diss, worst, axis=0)
            % Diss = np.delete(Diss, worst, axis=1)
            Memo = sum(Diss, 1);
            A(I(worst),:)=[];
        end
        temp=0;
    end
end