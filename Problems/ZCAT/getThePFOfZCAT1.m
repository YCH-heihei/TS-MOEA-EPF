% the sampled points 3-objective ZCAT  as follow:
% 适用于1-4，6-8
x=0:1/(100):1;
% x=0:1/(150):1;
[X,Y]=meshgrid(x);
% opt=gpuArray(PF1OfZCAT([X(:),Y(:)]));
%% 生成4万个点，接着反复删除与其它个体拥挤距离最小的解。
opt=PF1OfZCAT([X(:),Y(:)]);
problemsName='ZCAT2';
numOfPoints=size(opt,1);
M=size(opt,2);
dis=squareform(pdist(opt));
dis(logical(eye(size(opt,1))))=M;
[minDis,minIndex]=min(dis,[],2);

for i=1:10000
   % 随机生成一个新的点，将其加入到dis中，接着计算新的拥挤距离。
   newPoint=PF1OfZCAT([rand(1),rand(1)]);
   disOfNewPoints=pdist2(newPoint,opt);
   [minDisOfNewPoints,minIndexOfNewPoints]=min(disOfNewPoints);
   [minDisOfOpt,minIndexOfOpt]=min(minDis);
   if(minDisOfNewPoints>minDisOfOpt)
      opt(minIndexOfOpt,:)=newPoint;
%       try
            dis(minIndexOfOpt,[1:minIndexOfOpt-1,minIndexOfOpt+1:end])=disOfNewPoints([1:minIndexOfOpt-1,minIndexOfOpt+1:end]);
%       catch
          temp=0;
%       end
      dis(minIndexOfOpt,minIndexOfOpt)=M;
      dis(:,minIndexOfOpt)=dis(minIndexOfOpt,:);
%       disOfNewPoints(minIndexOfOpt)=M;
%       minIndexOfNewPoints=minIndexOfNewPoints+(minIndexOfNewPoints>minIndexOfOpt);
      if(minIndexOfNewPoints==minIndexOfOpt)
          [tempMinDis,tempIndex]=min(dis(minIndexOfOpt,:));
          minDis(minIndexOfOpt)=tempMinDis;
          minIndex(minIndexOfOpt)=minIndexOfOpt;
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
end
%% 将之保存到对应文件夹下。
folder = fullfile('.\Problems\Multi-objective optimization\ZCAT',problemsName);
[~,~]  = mkdir(folder);
% folder=folder(1:end-1);
file   = fullfile(folder,sprintf('%s_M%d_Num%d',problemsName,M,numOfPoints));
save([file,'.mat'],'opt');
temp=0;