% the sampled points 3-objective ZCAT  as follow:
x=0:1/(10):1;
% x=0:1/(150):1;
[X,Y]=meshgrid(x);
opt=gpuArray(PF1OfZCAT([X(:),Y(:)]));
%% 生成4万个点，接着反复删除与其它个体拥挤距离最小的解。
dis=gpuArray(squareform(pdist(opt)));
sortDis=gpuArray(zeros(size(opt,1)));
sortIndex=gpuArray(zeros(size(opt,1)));
for i=1:size(opt,1)
    [tempSortDisOfi,tempSortIndexOfi]=gpuArray(sort(dis(i,:),'ascend'));
    sortDis(i,:)=tempSortDisOfi;
    sortIndex(i,:)=tempSortIndexOfi;
end

problemsName='ZCAT1';
numOfPoints=101;
M=3;
while(size(opt,1)>numOfPoints)
    [~,p]=min(sortDis(:,2));
    q=sortIndex(p,2);
    % 删除最拥挤的个体。
    opt(q,:)=[];
    sortDis(q,:)=[];
    sortIndex(q,:)=[];
    tempSortIndex=gpuArray(zeros(size(opt,1),size(opt,1)));
    tempSortDis=gpuArray(zeros(size(opt,1),size(opt,1)));
    parfor i=1:size(opt,1)
       notSaveIndex=find(sortIndex(i,:)==q);
       tempSortDis(i,:)=sortDis(i,[1:notSaveIndex-1,notSaveIndex+1:end]);
       tempSortIndex(i,:)=sortIndex(i,[1:notSaveIndex-1,notSaveIndex+1:end]);
       tempSortIndex(i,tempSortIndex(i,:)>q)=tempSortIndex(i,tempSortIndex(i,:)>q)-1;
    end
    sortIndex=tempSortIndex;
    sortDis=tempSortDis;
end
%% 将之保存到对应文件夹下。
folder = fullfile('.\Problems\ZCAT',problemsName);
[~,~]  = mkdir(folder);
% folder=folder(1:end-1);
file   = fullfile(folder,sprintf('%s__M%d_Num%d',problemsName,M,Problem.M,numOfPoints));
save([file,'.mat'],'opt');
temp=0;