% the sampled points 3-objective ZCAT  as follow:
x=0:1/(100):1;

% x=0:1/(150):1;
[X,Y]=meshgrid(x);
% opt=gpuArray(PF1OfZCAT([X(:),Y(:)]));
PFName='PF16OfZCAT';
problemsName='ZCAT16';
%% 生成4万个点，接着反复删除与其它个体拥挤距离最小的解。
% opt=eval(strcat(PFName,'([X(:),Y(:)])'));
% opt=eval(strcat(PFName,'([X(:),Y(:)],3)'));
x=x';
% opt=eval(strcat(PFName,'(x)'));
opt=eval(strcat(PFName,'(x,3)'));
zMin=min(opt,[],1);
M=size(opt,2);
zMax=(1:M).^2;
range=zMax-zMin;


% 
% 

numOfPoints=10000;
[Fno,~]=NDSort(opt,1);
% PopDec=[X(:),Y(:)];
% notValidXAndY=PopDec(Fno~=1,:);
opt=opt(Fno==1,:);
normOfOpt=(opt-zMin)./(range);

dis=squareform(pdist(normOfOpt));
dis(logical(eye(size(opt,1))))=M;
[minDis,minIndex]=min(dis,[],2);
% 借助notValidXAndY;
% notValidXAndY=[];
for i=1:40000
   % 随机生成一个新的点，将其加入到dis中，接着计算新的拥挤距离。
%      newPoint=eval(strcat(PFName,'([rand(1),rand(1)])'));
%      inputDec=[rand(1),rand(1)];

     inputDec=[rand(1)];
     newPoint=eval(strcat(PFName,'(inputDec,3)'));
%      newPoint=eval(strcat(PFName,'(inputDec,3)'));
%      newPoint=eval(strcat(PFName,'(inputDec)'));
% %      newPoint=eval(strcat(PFName,'(inputDec)'));                                                                                         
%      newPoint=eval(strcat(PFName,'(inputDec,3)'));
     normOfNewPoint=(newPoint-zMin)./range;
     diff=normOfNewPoint-normOfOpt;

     temp=all(diff<0,2);
     dominateIndex=find(temp);
     notdominateIndex=find(~temp);
     % 存在被支配的个体，那么将被支配个体删掉。
     if(length(dominateIndex)~=0)
         opt(dominateIndex,:)=[];
         dis(dominateIndex,:)=[];
         dis(:,dominateIndex)=[];
         normOfOpt(dominateIndex,:)=[];
         diff(dominateIndex,:)=[];
         deleteStartIndex=min(dominateIndex);
         minIndex(dominateIndex)=[];
         minDis(dominateIndex)=[]; 

         % 更新与已删除个体最近的个体的最近距离。
         compared=repmat(dominateIndex,length(notdominateIndex),1); 
         temp=sum(minIndex==compared,2)>0;
         needUpdateIndex=find(temp);
         notNeedUpdateIndex=find(~temp);
         for j=1:length(needUpdateIndex)
             [tempDis,tempIndex]=min(dis(needUpdateIndex(j),:));
             minIndex(needUpdateIndex(j))=tempIndex;
             minDis(needUpdateIndex(j))=tempDis;
         end

         % 删除这个被支配个体后的新序号(减去删除的个体数)。
         needCorrectIndex=minIndex(notNeedUpdateIndex)>deleteStartIndex;
         minIndex(notNeedUpdateIndex(needCorrectIndex))=minIndex(notNeedUpdateIndex(needCorrectIndex))-length(dominateIndex);   
     end

     dominatedIndex=find(all(diff>0,2));
     % 新产生个体是被支配个体，那么直接不需要考虑。
     if(length(dominatedIndex)>0)
%         notValidXAndY=[notValidXAndY;inputDec];
        continue;
     end
     
     % 新个体被旧个体支配。
     disOfNewPoints=sum(diff.^2,2).^(1/2);
     if(size(opt,1)>=numOfPoints)
          [minDisOfNewPoints,minIndexOfNewPoints]=min(disOfNewPoints);
          [minDisOfOpt,minIndexOfOpt]=min(minDis);
          if(minDisOfNewPoints>minDisOfOpt)
              % 找到minDisOfNewPoints中的第二近节点，比较二者谁更小一些。
              p=minIndexOfOpt;
              q=minIndex(minIndexOfOpt);
              crow1=min(dis(p,[1:q-1,q+1:end]));
              crow2=min(dis(p,[1:p-1,p+1:end]));
              if(crow1>crow2)
                  minIndexOfOpt=q;
              else
                  minIndexOfOpt=p;
              end
              opt(minIndexOfOpt,:)=newPoint;
              normOfOpt(minIndexOfOpt,:)=normOfNewPoint;

              dis(minIndexOfOpt,[1:minIndexOfOpt-1,minIndexOfOpt+1:end])=disOfNewPoints([1:minIndexOfOpt-1,minIndexOfOpt+1:end]);

              dis(minIndexOfOpt,minIndexOfOpt)=M;
              dis(:,minIndexOfOpt)=dis(minIndexOfOpt,:);
      
              if(minIndexOfNewPoints==minIndexOfOpt)
                  [tempMinDis,tempIndex]=min(dis(minIndexOfOpt,:));
                  minDis(minIndexOfOpt)=tempMinDis;
                  minIndex(minIndexOfOpt)=tempIndex;
              else
                  minDis(minIndexOfOpt)=minDisOfNewPoints;
                  minIndex(minIndexOfOpt)=minIndexOfNewPoints;
              end
          
              % 找到minIndex为minIndexOfopt的个体
              needUpdatedIndex=minIndex==minIndexOfOpt;
              [newMinDisOfNeeUpdated,newTempIndexOfNeeUpdated]=min(dis(needUpdatedIndex,:),[],2); 
              minDis(needUpdatedIndex)=newMinDisOfNeeUpdated;
              minIndex(needUpdatedIndex)=newTempIndexOfNeeUpdated;
          
              % 更新minDis
              needUpdatedIndex=(minDis>dis(minIndexOfOpt,:)');
              
        %       [newMinDisOfNeeUpdated,newTempIndexOfNeeUpdated]=min(dis(needUpdatedIndex,:),minDis()); 
              minDis(needUpdatedIndex)=dis(needUpdatedIndex,minIndexOfOpt);
              minIndex(needUpdatedIndex)=minIndexOfOpt;
          end
     else
          % 将新个体直接追加到opt中。
          opt(end+1,:)=newPoint;
          normOfOpt(end+1,:)=normOfNewPoint;
%           try
             dis(:,end+1)=disOfNewPoints;
%           catch
%              temp=0;
%           end
            disOfNewPoints(end+1)=M;
          dis(end+1,:)=disOfNewPoints;
          [tempDis,tempIndex]=min(disOfNewPoints);
          minIndex(end+1)=tempIndex;
          minDis(end+1)=tempDis;
          
          
          % 判断新加入节点是否改变了种群节点最近距离，若改变了则将其最近节点进行更新。
          needUpdateIndex=find(minDis(1:end-1)<disOfNewPoints(1:end-1));
          minIndex(needUpdateIndex)=size(opt,1);
          minDis(needUpdateIndex)=disOfNewPoints(needUpdateIndex);
       
     end
end

if(M==2)
    [~,orderIndex]=sortrows(opt);
    opt=opt(orderIndex,:);
end
%% 将之保存到对应文件夹下。
folder = fullfile('.\Problems\Multi-objective optimization\ZCAT',problemsName);
[~,~]  = mkdir(folder);
% folder=folder(1:end-1);
file   = fullfile(folder,sprintf('%s_M%d_Num%d',problemsName,M,numOfPoints));
% save([file,'.mat'],'opt','notValidXAndY');
save([file,'.mat'],'opt');
temp=0;