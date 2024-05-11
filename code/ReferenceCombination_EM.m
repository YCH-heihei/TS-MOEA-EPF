function [Ruq,net] = ReferenceCombination_EM(Ru,net)
% function [Ruq,net,Ru] = ReferenceCombination_EM(Ru,net,EM)
% Combine uniform reference vectors and nodes in GNG (Algorihtm 4)
% Ru: 初始均匀分布的参考向量;
% net: GNG;
% GNG: 具有自动划分聚类以及记录形状的特点，能否与LMPFG中的模型进行结合，从而适应更复杂的前沿？ 

%--------------------------------------------------------------------------
% Copyright Yiping Liu
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------

    numNode = size(net.NodeS,1);
    % NodeP: % Node mapped to hyperplane 
    % NodeS: Expanded node 
    net.NodeP = net.NodeS; 
%     M=size(ExtremeVector,2);
%     if size(net.NodeS,1)<=1
%         Ruq=Ru;
%          ratio=sum(ExtremeVector(1:M,:),2);
%          expandRatio=1.^-2*(repmat(ratio,1,M));
%          ExtremeVector(1:M,:)=ExtremeVector(1:M,:)+eye(M).*expandRatio;
%          tempExtremeVector=ExtremeVector(1:M,:)-ExtremeVector(1:M,:).*logical(eye(M));
%          tempExtremeVector(tempExtremeVector<10.^(-6))=10.^-6;
%          ExtremeVector(1:M,:)=ExtremeVector(1:M,:)-(expandRatio).*(tempExtremeVector./repmat(sum(tempExtremeVector,2),1,M));
%          ExtremeVector(ExtremeVector<10.^-6)=10.^-6;
%          vector=ExtremeVector./(repmat(sum(ExtremeVector,2),1,M));
%          ExtremeVector(ExtremeVector<10.^(-6))=10.^(-6);
%          EM=vector;
% %          q=0.00;
% %          [~,emIndex]=max(net.NodeP,[],1);
%         return
%     end
    
    
    %% Map nodes to hyperplane(各维度相加为1)，将其投影到平面为1的平面上。
    for i = 1:numNode
       x = 1./sum(net.NodeS(i,:));
       net.NodeP(i,:) = net.NodeS(i,:)*x; 
    end    
    
      %% 除此之外，还需要对极值向量让其有向外扩的趋势，这里采用的方法是让最大的目标上+10^-2。
     %% 而其余目标减去其10^-2*其余目标原目标的占比，这样对应参考更倾向于选择当前极值点对应的方向。
     
%      ratio=sum(ExtremeVector(1:M,:),2);
%      expandRatio=1*10.^-4*(repmat(ratio,1,M));
%      ExtremeVector(1:M,:)=ExtremeVector(1:M,:)+eye(M).*expandRatio;
%      tempExtremeVector=ExtremeVector(1:M,:)-ExtremeVector(1:M,:).*logical(eye(M));
%      tempExtremeVector(tempExtremeVector<10.^(-6))=10.^(-6);
%      ExtremeVector(1:M,:)=ExtremeVector(1:M,:)-(expandRatio).*(tempExtremeVector./repmat(sum(tempExtremeVector,2),1,M));
%      ExtremeVector(ExtremeVector<10.^(-6))=10.^(-6);
%      vector=ExtremeVector./(repmat(sum(ExtremeVector,2),1,M));
%      ExtremeVector(ExtremeVector<10.^(-6))=10.^(-6);
%      q=0.01;
%      [cosValue,emIndex]=max(1-pdist2(vector,net.NodeP,'cosine'),[],2);
%      EM=zeros(length(emIndex),M);
%      for i=1:length(emIndex)
%         [~,index]=max(net.NodeP,[],1);
%         [cosValue,index]=max(1-pdist2(vector(i,:),net.NodeP,'cosine'));
%         方案1:
%         net.Node(emIndex(i),:)=q*net.Node(emIndex(i),:)+(1-q)*ExtremeVector(i,:);
%         net.NodeP(emIndex(i),:)=q*net.Node(emIndex(i),:)./(sum(net.Node(emIndex(i),:)));
%         方案2: 
%         net.NodeP(emIndex(i),:)=net.Node(emIndex(i),:)./sum(net.Node(emIndex(i),:));
        %% ./(sum(net.NodeP(index,:)));
%         EM(i,:)=q*net.NodeP(emIndex(i),:)+(1-q)*vector(i,:);
%         EM(i,:)=EM(i,:)./sum(EM(i,:),2);
%         temp=0;
%      end
%     EM=net.NodeP(emIndex,:);
%     temp=0;
    
    %% Average distance ammong nodes
    % 对应两组节点若存在边，那么此边的权重为此两点之间的距离(问题就在于这里，即使某两个相邻点很远，并且很久之前存在点，但是这里误以为没有点，从而认为这个区域不存在个体，需要不需要引入参考向量)
    % 需要针对我们的问题来修改阈值，因为我们使用的方法可能会将相距稍远的个体相连接。
    % 原来的目的在于通过将累积误差积累到边上而记录数据点出现在此处。
    % 但没想到会因为边平均距离的增加而导致阈值增大，从而即使存在空虚的区域也不会引入参考向量。
    % 但是又不能引入过多的外部参考向量，例如本来有些区域就是被其它方向支配的，但还是引入了这个方向。
    % 怎么办，即需要保留探索GNG节点以外方向个体的优化,又要避免在无效方向上关联个体。
    % 如果把阈值设的太小，那么会引入过多初始化的方向，如果为这些无效方向关联一个个体，而这个方向又不能收敛到针对拍雷托前沿上。


%   
     %% Choose the smaller one
%      AvgDis = min(AvgDis,MinDis);
% 
%      %% Remove some reference vectors in Ru which are too close to nodes
          
     tempChoose=zeros(1,size(Ru,1));
     if(size(net.NodeP,1)>=1 && size(net.NodeP,1)<net.maxNode)
        Distance3 = pdist2(Ru,net.NodeP,'cosine');
        %% 限制引入的外部节点数的数目。
        Distance3=min(Distance3,[],2);
%         Distance3=sort(Distance3,2,'ascend');
        [~,order]=sort(Distance3,'descend');
%         [~,order]=sort(sum(Distance3(:,min(1:min(size(Ru,2),size(Distance3,2)))),2),'descend');
        tempChoose(order(1:size(Ru,1)-size(net.NodeP,1)))=1;
        Choose=tempChoose==1;
         Ruq=Ru(Choose,:);
     elseif(size(net.NodeP,1)==0) 
         Ruq=Ru;
     else 
         Ruq=[];
     end
     temp=0;
    
     
     
   
end