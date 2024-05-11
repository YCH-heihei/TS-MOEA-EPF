function Choose = ArchiveUpdate_MD(Data,N1,net,Zmin,zMax,EM)
%% 极值点方案2 function [Data,nND] = ArchiveUpdate_MD(Data,N1,N2,Z1,Z2,Zmin,EM)
% Update Input Signal Archive in DEA-GNG
% Data: 上一次的用于训练的信号，以及当前代的目标向量;
% N: 训练数据的上限;
% Z1: 用于选择的参考向量;
% Z2: GNG中展开后的节点;
% Zmin: 理想点;

% Data: 

%--------------------------------------------------------------------------
% Copyright Yiping Liu
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------
    
    %% Delete duplicate
    
%     [N,M]=size(Data);
    %% 删除那些没有显著优势但是值却很差的点(归一化后，加上细微噪声就被其它个体支配掉的解)。

    
    %% Non-dominated sorting
%     FrontNo = NDSort(Data,N1);
    % 优先选择第一层，因为第一层前沿最接近真实前沿，为什么第一层最接近真实前沿?
    % 既然只需要选择第一层那么可以直接判断个体是否与其它个体相比具有显著非支配关系。

     % 不足的话，直接为第一层的个体。

    nND = size(Data,1);
    Choose=true(1,nND);
%     A=sum(Data.^(2),2);
%     B=mean(A,2);
%     Data(A>10*B,:)=[];
    %% Selection
    if nND > N1
        % 对支配等级为1的个体进一步进行筛选。
        % 极值点方案2: Choose = LastSelection(Data,N1,Z1,Z2,Zmin,EM);
        Choose = LastSelection(Data,N1,net,Zmin,EM);
        Data = Data(Choose',:); 
        nND = N1;      
    end
end

%% 这个算法所选择的个体全都为同一层的个体, 然后针对参考向量均分保留个体,而非保留在一处的个体.
%% 从而期待对GNG中的所有节点均有训练,避免某些节点长期未更新而被删去,这是一方面.
%% 另外一方面是避免某些地方更新得很好,其它方面更新不是很好.
%% 这个算法所选择的个体全都为同一层的个体, 然后针对参考向量均分保留个体,而非保留在一处的个体.
%% 从而期待对GNG中的所有节点均有训练,避免某些节点长期未更新而被删去,这是一方面.
%% 另外一方面是避免某些地方更新得很好,其它方面更新不是很好.
function Choose = LastSelection(PopObj,K,net,Zmin,EM)
%% 极值点方案2: function Choose = LastSelection(PopObj,K,Z1,Z2,Zmin,EM)
% Select part of the solutions in the last front
     % 使用聚类算法，划分为K个聚类，然后选择那些与其它个体差异较大的。
%         [Zmax,maxIndex] = max(PopObj,[],1);
        [Zmax,maxIndexOfPop] = max(PopObj,[],1); 
%         [~,minIndexOfPop]= min(PopObj,[],1);
%         [Zmax,~] = max(PopObj,[],1);
        [N,M]=size(PopObj);
        PopObj=(PopObj - repmat(Zmin,N,1))./repmat(Zmax-Zmin,N,1);  


%     % 利用关联的参考向量中垂直距离最小的来找到极值点。
%     %% 不能使用上面的切比雪夫(只有当对应极值与其它目标差异非常大时，使用切比雪夫才有效)来找到极值点，修改为使用垂直距离的方式  

    [cosValue,pi]=max(1-pdist2(PopObj,EM,'cosine'),[],2);
    normOfObjs=sum(PopObj.^2,2).^(1/2);
    VerticalDis=normOfObjs.*((1-cosValue.^(2)).^(1/2));
    pjDis=normOfObjs.*((cosValue.^(2)).^(1/2));
    
    maxIndex=[];
    for i=1:M
         current=find(pi==i);
         if(length(current)~=0)
            [~,index]=min(VerticalDis(current)+0.05*pjDis(current));
            maxIndex(end+1)=current(index);
         end
    end
    
% 方案3: 求出极值点构成的超平面。
      
%%  聚类算法的方式。
%     [~,maxNormIndex]=max(normOfObjs);
%     saveIndex=maxNormIndex;
%       saveIndex=maxIndexOfPop;
%     dis=pdist_cosineAndEuclidean(PopObj(1:K,:));
    % 先选择前K个，然后反复从后N-K加入后立马进行筛选。
%     dis=squareform(pdist(PopObj(1:K,:)));
      allDis=squareform(pdist(PopObj));
      allDis(logical(eye(N)))=M;
      cosineDis=squareform(pdist(PopObj,'cosine'));
      cosineDis(logical(eye(N)))=M;
%     allDis=squareform(pdist(PopObj)).*(squareform(pdist(PopObj,'cosine')));% 计算所有个体的距离。 
%     allDis=(squareform(pdist(PopObj,'cosine')));
       
    
%     minDisFromAllDis=min(allDis,[],2);
%     [~,orderOfMDFAD]=sort(minDisFromAllDis,'descend');
%     Choose  = false(1,N);
    % 起点设为与原点距离最大的点。
    
%     saveIndex=maxIndex; % 新方案，先选择极值点，再反复选择与已选择节点距离最近距离最小的哪一个。
%     notSaveIndex=setdiff(1:N,saveIndex);
%     closestDis= min(allDis(notSaveIndex,saveIndex),[],2); % 与其它个体的最近距离。
%     while(length(saveIndex)<K)
% %         notSaveIndex=setdiff(1:N,saveIndex);
% %         minDis=min(allDis(notSaveIndex,saveIndex),[],2);
%         [~,appendIndex]=max(closestDis);
%         saveIndex=[saveIndex,notSaveIndex(appendIndex)];
%         notSaveIndex(appendIndex)=[];
%         closestDis(appendIndex)=[];
% %         closestTosaveIndex=[];
%         newDis=allDis(notSaveIndex,saveIndex(end));
%         changeIndex=newDis<closestDis;
%         closestDis(changeIndex)=newDis(changeIndex);
%     end
    
    saveIndex=maxIndex; % 新方案，先选择极值点，再反复选择与已选择节点距离最近距离最小的哪一个。
    notSaveIndex=setdiff(1:N,saveIndex);
%     cosineDis=squareform(pdist(PopObj,'cosine'));
%     cosineDis(logical(eye(N)))=M;
%     closestDis= min(allDis(notSaveIndex,saveIndex),[],2); % 与其它个体的最近距离。
    closestDis= min(cosineDis(notSaveIndex,saveIndex),[],2); % 与其它个体的最近距离。
    while(length(saveIndex)<K)
%         notSaveIndex=setdiff(1:N,saveIndex);
%         minDis=min(allDis(notSaveIndex,saveIndex),[],2);
        [~,appendIndex]=max(closestDis);
        saveIndex=[saveIndex,notSaveIndex(appendIndex)];
        notSaveIndex(appendIndex)=[];
        closestDis(appendIndex)=[];
%         closestTosaveIndex=[];
%         newDis=allDis(notSaveIndex,saveIndex(end));
        newDis=cosineDis(notSaveIndex,saveIndex(end));
        changeIndex=newDis<closestDis;
        closestDis(changeIndex)=newDis(changeIndex);
    end
    saveIndex=[saveIndex,1];
    dis=[];

    saveIndex=[saveIndex,1];
    dis=[];
    cosineDisOfSave=[];
    for i=1:K
        index=[saveIndex(1:i-1),saveIndex(i),saveIndex(i+1:K)];
        dis=[dis;allDis(saveIndex(i),index)];
        cosineDisOfSave=[cosineDisOfSave;cosineDis(saveIndex(i),index)];
    end 
    I=eye(K+1);
    for i=1:K
       for j=i+1:K
            temp=PopObj(saveIndex(j),:)-PopObj(saveIndex(i),:); 
            I(i,j)=(max(temp));
            I(j,i)=(max(-temp));
       end
    end
%     saveIndex=1:K+1; % 原有方案。
%     saveIndex=orderOfMDFAD(1:K+1); 
    
    insertIndex=K+1;
    nodeDis=squareform(pdist(net.NodeS));
%     nodeDis=pdist_cosineAndEuclidean(net.NodeS);
    maxDis=1*sum(nodeDis(:).*net.edge(:))./(sum(net.edge(:)));
    [minDis,minIndex]=min(cosineDisOfSave,[],2);
%     [minDis,minIndex]=min(dis,[],2);
%     cosineDis
    for i=notSaveIndex
        saveIndex(insertIndex)=i;
        notIndexIndex=[1:insertIndex-1,insertIndex+1:K+1];
        otherIndex=saveIndex(notIndexIndex);
        % 与现有其它个体的欧式距离。
        % 将新个体与现有个体的距离更新到对应的距离矩阵中。
%         tempDis=pdist2(PopObj(otherIndex,:),PopObj(i,:)).*(pdist2(PopObj(otherIndex,:),PopObj(i,:),'cosine').^(1/2));
        % 列赋值。
%         tempDis=pdist2(PopObj(otherIndex,:),PopObj(i,:));
        tempDis=allDis(i,otherIndex); % 插入节点与现有节点的距离。
        index=[1:insertIndex-1,insertIndex+1:K+1];
        dis(index,insertIndex)=tempDis;
        % 列赋值。
        dis(insertIndex,index)=tempDis;
        dis(insertIndex,insertIndex)=M;
        tempDis=cosineDis(i,otherIndex);
        cosineDisOfSave(index,insertIndex)=tempDis;
        cosineDisOfSave(insertIndex,index)=tempDis;
        cosineDisOfSave(insertIndex,insertIndex)=M;
        temp=PopObj(otherIndex,:)-PopObj(i,:);

        I(insertIndex,index)=max(temp,[],2);
        I(index,insertIndex)=max(-temp,[],2);
        I(insertIndex,insertIndex)=1;
        
        % 只需要找到M个大于最小值的即可，不需要更新排序关系了，只需要记录个体与其它哪个个体的距离最小。在找到大于。
       
        % 清空。
        
         [minDisOfInsert,minDisOfInsertIndex]=min(cosineDisOfSave(insertIndex,:));
%          [minDisOfInsert,minDisOfInsertIndex]=min(dis(insertIndex,:));
         minDis(insertIndex)=minDisOfInsert;
         minIndex(insertIndex)=minDisOfInsertIndex;
         % 因insertIndex被删除，而必要重新跟新的结点。
         needUpdateIndex=find(minIndex==insertIndex);
         [minDisOfUpdate,minDisOfUpdateIndex]=min(cosineDisOfSave(needUpdateIndex,:),[],2);
%          [minDisOfUpdate,minDisOfUpdateIndex]=min(dis(needUpdateIndex,:),[],2);
         % changeIndex: 需要重新的结点。
         minDis(needUpdateIndex)=minDisOfUpdate;
         minIndex(needUpdateIndex)=minDisOfUpdateIndex;
         
%          try
            changeIndex=find(minDis>cosineDisOfSave(:,insertIndex));
%             changeIndex=find(minDis>dis(:,insertIndex));
%          catch
            temp=0;
%          end
         minDis(changeIndex)=cosineDisOfSave(changeIndex,insertIndex);
%          minDis(changeIndex)=dis(changeIndex,insertIndex);
         minIndex(changeIndex)=insertIndex;
         
         [~,rb]=min(minDis);
%         maxDis=max(0.02*max(max(dis)),2*mean(minDis));

        
         ra=minIndex(rb);
%         try
%             indexOfRa=dis(ra,:)<maxDis;
%             indexOfRb=dis(rb,:)<maxDis;
         indexOfRa=find(dis(ra,:)<maxDis);
         indexOfRb=find(dis(rb,:)<maxDis);
%         catch
            
         temp=0;
        [sortDisOfRa,sortIndexOfRa]=sort(dis(ra,indexOfRa),'ascend');
        [sortDisOfRb,sortIndexOfRb]=sort(dis(rb,indexOfRb),'ascend');
        num=1*M;
        if(length(sortDisOfRa)>length(sortDisOfRb))
            num=length(sortDisOfRb);
        else
            num=length(sortDisOfRa);
        end
        num=min(num,1*M);

%         end 
%         IOfRa=sum(I(ra,indexOfRa).^(2).*dis(ra,indexOfRa).^(1));
%         IOfRb=sum(I(rb,indexOfRb).^(2).*dis(rb,indexOfRb).^(1));
        IOfRa=sum(I(ra,indexOfRa(sortIndexOfRa(1:num))).^(1).*sortDisOfRa(1:num).^(1));
        IOfRb=sum(I(rb,indexOfRb(sortIndexOfRb(1:num))).^(1).*sortDisOfRb(1:num).^(1));

        if(IOfRa>IOfRb)
            insertIndex=rb;
        else
            insertIndex=ra;
        end
%         if(i==N-M)
%             
%                 temp=0;
%         end
        temp=0;
    end 
    %% 旧方案。
%     saveIndex(insertIndex)=maxIndex(1); % 用必选节点来清空删除的节点。
    saveIndex(insertIndex)=[];
    Choose(saveIndex)=true;
    Choose(maxIndex)=true;
%     Choose(minIndexOfPop)=true;
    %% 旧方案。   
end