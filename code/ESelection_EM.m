function [Population,FrontNo,net,Zmax] = ESelection_EM(Population,FrontNo,MaxFNo,N,M,Ru,net,Zmin,Zmax,theta,EM)
% function [Population,FrontNo,crd,EM] = ESelection_EM(Population,N,Ruq,Rnode,theta,Zmin,EM)
% The environmental selection of DEA-GNG
% Population: 待筛选的种群。
% N: 待留下的个体数。
% Ruq: 参考向量。
% Rnode: 当前GNG中的节点。
% thete: PBI的权重。
% Zmin: 理想点。


%--------------------------------------------------------------------------
% Copyright Yiping Liu
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------

    %% Non-dominated sorting
%     [numOfIndividuals,M]=size(Population.objs);
%        FrontNo=[FrontNo,ones(1,length(AS))];
%        Population=[Population,AS];
%        [~,ia,~]=unique(Population.objs,'rows');
%        Population= Population(ia);
%        FrontNo= FrontNo(ia);
       MaxFNo=1;
       while(sum(FrontNo<=MaxFNo)<N)
          MaxFNo=MaxFNo+1;
       end
%     temp=0;
    % 优先选择除倒数第一层的前沿。
%         F1 = Population(FrontNo==1);
        
%         sumTemp=sum(temp,2);
%         maxIndex=[];
%         for i=1:M
%             tail=1*prod(temp(:,[1:i-1,i+1:M]).*2,2)/M;
%             scoreOfMax=-temp(:,i)./sumTemp+tail;
%             [~,minIndex]=min(scoreOfMax);
%             maxIndex=[maxIndex,minIndex];
%         end
        
        
%       cosineToLine=pdist2(temp,ones(1,M),'cosine');
%       [~,maxIndex]=max(cosineToLine);
%       cosineToLine(maxIndex)=[];
%       notInMaxIndex=setdiff(1:length(temp),maxIndex);
%       minCosDis=min(pdist2(temp(notInMaxIndex,:),temp(maxIndex(end),:),'cosine'),cosineToLine);
% %       minCosDis=pdist2(temp(notInMaxIndex,:),temp(maxIndex(end),:),'cosine')+cosineToLine;
%       while length(maxIndex)<M
%           [~,appendIndex]=max(minCosDis);
%           maxIndex=[maxIndex,notInMaxIndex(appendIndex)];
%           minCosDis(appendIndex)=[];
%           notInMaxIndex(appendIndex)=[];
%           minCosDis=min(minCosDis,pdist2(temp(notInMaxIndex,:),temp(maxIndex(end),:),'cosine'));
%       end
%         ExtremeVectors=temp(maxIndex,:);
        Next = FrontNo < MaxFNo;
        Last   = find(FrontNo==MaxFNo);
        [Ruq,net] = ReferenceCombination_EM(Ru,net);    
        
        
        Choose = LastSelection(Population(Next).objs,Population(Last).objs,N-sum(Next),Ruq,net.NodeP,theta,Zmin,Zmax,EM);
    
        %% Select the solutions in the last front
        Next(Last(Choose)) = true;

        %% Population for next generation
        Population = Population(Next);
        FrontNo    = FrontNo(Next); 
        temp=0;
end
function Choose = LastSelection(PopObj1,PopObj2,K,Ruq,Rnode,theta,Zmin,Zmax,EM)
%%  方案1 :function [Choose, crd] = LastSelection(PopObj1,PopObj2,K,Ruq,Rnode,theta,Zmin,Zmax,EM)
% Select part of the solutions in the last front
% PopObj1: 已选择的个体。
% PopObj2: 待选择的个体。
% K: 需选择的个体数。
% Ruq: 参考向量。
% Rnode: GNG中的节点。
% theta: PBI权重。
% Zmin: 理想点。
% Zmax: 极差点。


    %% Initialize
    PopObj = [PopObj1;PopObj2]; %Candidate solutions
    R = [Ruq;Rnode]; % Reference Vectors   
    
    cosineOfR = 1 - pdist2(R,R,'cosine');
    
%     cosineOfR = 
    cosineOfR (logical(eye(length(cosineOfR)))) = 0;
%     [sortCos,~]=sort(cosineOfR,2);
    angle = acos(cosineOfR);
    angle=sort(angle,2,'ascend');
%     gamma  = min(acos(cosineOfR),[],2);
    
%     gamma  = min(acos(cosineOfR),[],2); % 参考向量与其它参考向量的最小夹角，用以作为APD值标准化所使用的值。
    
    
    [N,M]  = size(PopObj);
%     gamma  = min(angle(:,M),[],2);
    gamma  = mean(angle(:,1:M),2);
    N1     = size(PopObj1,1);
    N2     = size(PopObj2,1);
    NR     = size(R,1);
    NR1    = size(Ruq,1);   
        
    %% Normalization [0-1]
    PopObj = (PopObj - repmat(Zmin,N,1))./repmat(Zmax-Zmin,N,1);
   
    %% Scalarizing function value
    %% Associate each solution with one reference vector
    % Calculate the distance of each solution to each reference vector
    Cosine= 1 - pdist2(PopObj,R,'cosine');
    Angle= acos(Cosine);
%     NormP = sqrt(sum((PopObj.*(Zmax-Zmin)).^2,2));  %
%     APD中使用的，虽然在MaF7上表现上很好，但感觉不合理。
    scaled=mean(PopObj,1);
    scaled=(scaled./min(scaled)).^(2);
%     scaled=ceil();
    NormP = sqrt(sum((PopObj.*scaled).^2,2));
%       NormP = sqrt(sum((2*PopObj).^2,2));
%     NormP=2*min(PopObj,[],2)+max(PopObj,[],2);
    [d2,pi] = min(Angle,[],2);

    %% Select solutions
    Choose  = false(1,N2);
    %% Calculate the number of associated solutions except for the last front of each reference point
    % 以扰动来选择极值点。
    if(N1==0)
        [cosValue,piOfEM]=max(1-pdist2(PopObj,EM,'cosine'),[],2);
        normOfObjs=sum(PopObj.^2,2).^(1/2);
        VerticalDis=normOfObjs.*((1-cosValue.^(2)).^(1/2));
        pjDis=normOfObjs.*((cosValue.^(2)).^(1/2));
        maxIndex=[];
        for i=1:M
             current=find(piOfEM==i);
             if(length(current)~=0)
                % 40倍的的投影距离。
                [~,index]=min(VerticalDis(current)+0.025*pjDis(current));
                maxIndex(end+1)=current(index);  
             end
        end
% %       % 方案2: 选择极轴方向的极值点。
%         Choose(maxIndex)=1;
    end
       
    %% Calculate scalarizing functions (PBI)
    orderOfPop=ones(1,N2)*M;
%     orderOfPop(1:N1)=0;
    g = zeros(1,N2); 
    for i=1:NR
        % 找到当前参考向量关联的个体未选择曾关联的个体。
        current=find(pi(N1+1:N)==i);
        g(current)=NormP(current+N1).*(1+1*theta*(d2(current+N1))./(gamma(i)));
        
        [~,currentOrder]=sort(g(current),'ascend');
        % 还需要反序。
%         tempIndex=1:length(current);
        orderOfPop(current(currentOrder))=1:length(current);
        
    end 
    if(N1==0)
        g(maxIndex)=0;
    end
    for i=1:NR
        % 找到当前参考向量关联的个体未选择曾关联的个体。
        current=find(pi(N1+1:N)==i);
%         g(current)=NormP(current+N1).*(1+1*theta*(d2(current+N1))./(gamma(i)));
        
        [~,currentOrder]=sort(g(current),'ascend');
        % 还需要反序。
%         tempIndex=1:length(current);
        orderOfPop(current(currentOrder))=1:length(current);
    end 
    
%     if(sum(orderOfPop<=1)>K+N1)
%         % 待选择个体将从第一层中进行(N1为0)
%         orderOfPop(maxIndex)=1;
%         %% 反复选择与已经选择个体最小距离最大的那一个加入到已选择个体当中。
%         current=find(orderOfPop<=1);
%         allDis=squareform(pdist(PopObj));
%         allDis(logical(eye(N)))=M;
%         saveIndex=maxIndex; % 新方案，先选择极值点，再反复选择与已选择节点距离最近距离最小的哪一个。
%         notSaveIndex=setdiff(current,maxIndex);
%         closestDis= min(allDis(notSaveIndex,saveIndex),[],2);
%         while(length(saveIndex)<K+N1)
%                 [~,appendIndex]=max(closestDis);
%                 saveIndex=[saveIndex,notSaveIndex(appendIndex)];
%                 notSaveIndex(appendIndex)=[];
%                 closestDis(appendIndex)=[];
%         %         closestTosaveIndex=[];
%                 newDis=allDis(notSaveIndex,saveIndex(end));
%                 changeIndex=newDis<closestDis;
%                 closestDis(changeIndex)=newDis(changeIndex);
%         end
%         Choose(saveIndex)=true;
%     else

%         if(N1==0)
%             for i=maxIndex
%                 C=find(pi==pi(i));
%                 orderOfPop(C)=orderOfPop(C)+1;
%                 orderOfPop(i)=0;
%             end
%         end
        [~,order]=sort(orderOfPop,'ascend');
        Choose(order(1:K))=true;
%     end 

   
    
end