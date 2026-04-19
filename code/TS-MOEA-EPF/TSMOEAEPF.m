classdef TSMOEAEPF < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% Decomposition based evolutionary algorithm guided by growing neural gas
% aph ---   0.1  --- Parameter alpha
% eps --- 0.314 --- Parameter epsilon


%------------------------------- Reference --------------------------------
% Y. Liu, H. Ishibuchi, N. Masuyama, and Y. Nojima, Adapting reference
% vectors and scalarizing functions by growing neural gas to handle
% irregular Pareto fronts. IEEE Transactions on Evolutionary Computation,
% 2020, 24(3): 439-453.
%--------------------------------------------------------------------------
% Copyright Yiping Liu
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            tic;
            %% Parameter setting
            [aph,eps] = Algorithm.ParameterSet(0.1,0.314); 
            
            %% DEA Initialization
            [Ru,Problem.N] = UniformPoint(Problem.N,Problem.M);	% Uniform reference vectors
            Population     = Problem.Initialization();          % Random population
            zMin           = min(Population.objs,[],1);         % Ideal Point  
            AS             = [];                                % Input Signal Archive
            Ruq            = Ru;                                % Refernce vectors in Ru for selection
            PopObjs=Population.objs;
            [FrontNo,~]    = NDSort(Population.objs,Problem.N); % First Fitness for the first mating selection
            firstIndex=find(FrontNo==1);
            PopObjs=PopObjs(firstIndex,:);
            
            zMax=max(PopObjs,[],1);
            
            % 最小点的存在对离群点的判断有一定影响，所以每代还是强制保留为好，至。
            [~,minPointIndex]=min(sum((PopObjs-zMin)./(zMax-zMin),2));
            minOfPoints=PopObjs(minPointIndex,:);
            % 存档的初始化为非被支配个体。
            AS=Population(firstIndex);
            M=Problem.M;
%             AS=[PopObjs];
%             zMax=max(PopObjs,[],1);
            crd            = zeros(1,Problem.N);                % second Fitness for the first mating selection   
            MaxGen         = ceil(Problem.maxFE/Problem.N);    	% Maximum Generation
%             FrontNo = NDSort(Population.objs,Problem.N);
            
            %% GNG Initialization
%             spand=floor(1/20*Problem.M*Problem.N);
%             spand=ceil(1/4*Problem.N);      
%             ArchiveSize1 = ceil(1/2*Problem.M)*Problem.N;
%             ArchiveSize2 = ceil(1/2*Problem.M)*Problem.N;
            ArchiveSize1 = ceil(2)*ceil(Problem.N+3);
            ArchiveSize2 = ceil(2)*ceil(Problem.N+3);
            NoG = ceil(aph*MaxGen);                   % Number of generations of Not Training GNGNumber of iterations to train GNG per Generation(针对每轮训练需要对数据进行多少次重复训练)  
            GNGnet.maxIter = 1;                 % Number of iterations to train GNG per Generation  
%             GNGnet.maxAge = (Problem.N);          % Maximum cluster age 
            GNGnet.maxAge = (6*Problem.M); 
%             GNGnet.maxAge = (Problem.M).^2;
            GNGnet.maxNode = ceil(1*Problem.N);         % Max number of nodes
            GNGnet.lambda = 0.05*Problem.N;      % Cycle for topology reconstruction   
            GNGnet.hp = [];                     % Hit point of node 
            GNGnet.maxHP = 2*ArchiveSize2;       % Max HP of node 
            GNGnet.Node = [];                   % Node
            GNGnet.NodeS = [];                  % Expanded node 
            GNGnet.NodeP = [];                  % Node mapped to hyperplane 
            GNGnet.Err = zeros(2,2);            % Error(修改，将误差累积到边，而不是节点上) 
            GNGnet.edge = zeros(2,2);           % Edge between nodes 
            GNGnet.age = zeros(2,2);            % Age of edge 
            GNGnet.epsilon_a = 0.2;             % Learning coefficient
            GNGnet.epsilon_nb = 0.02;           % Learning coefficient of neighbor
            GNGnet.alpha = 0.5;                 % Nodes r1max and r2max error reduction constant
            GNGnet.delta = 1/5.^(1/ArchiveSize1);                 % Error reduction coefficient  
            trainedFlag=0;
            for i=1:20
                GNGnet = GNGUpdateByDirection(AS.objs,GNGnet,Problem.M+1,zMin,eye(Problem.M));
            end

            axisDiag=ones(1,M)./(M);
            Vector1=zeros(1,M);
            Vector1(1)=1;
            alphaAngle=acos(1-pdist2( Vector1, axisDiag,'cosine'));
            convexDgree=max(1/((0.125*2*M)/2),4/(1.4/1.1));
            concaveDgree=min((0.875*2*M/2),1.8./(1.4/2));
            
            EM=eye(Problem.M);
            rate=ones(1,Problem.M)/M;
            %% Optimization
           
%             [outlierIndex,minOfPoints] = determinedOutlier(AS.objs,AS.objs,zMin,);
            while Algorithm.NotTerminated(Population)


              %% 原始方案。
                                          
              
%               % GNG-based adaptation(当前代属于可训练GNG的代)
                if ceil(Problem.FE/Problem.N) <= MaxGen - NoG
                    %% 原始方案。
%                     try
                        MatingPool=TournamentSelection(2,Problem.N,FrontNo,crd); % crd: 反映拥挤度，还是挺有用的，尤其是在边界处的点。
%                      MatingPool=TournamentSelection(2,Problem.N,FrontNo);
%                     catch
%                         temp=0;
%                     end
%                     [~,idealIndex]=min(crd);
%                     if(sum(MatingPool==idealIndex)==0)
%                          
%                          MatingPool(1)=idealIndex;
%                     end
                   
                    FirstIndex=find(FrontNo==1);
                    PopObjs=Population.objs;
                    refPop=PopObjs(FirstIndex,:);
%                     refPop=PopObjs(FrontNo==min(FrontNo),:);
                    objs=PopObjs(FirstIndex,:);
                    [zMax,maxIndex]=max(objs,[],1);

                    Offspring  = OperatorGA(Problem,Population(MatingPool));
                     if(Problem.FE>=16340)
                            % rng(20000);
                            temp=0;
                     end
                     if(Problem.FE>=209)
                            temp=0;
                     end
                    zMin= min([zMin;Offspring.objs],[],1);
                    % Input Signal Archive Update
                    tempPopulation=[Offspring,Population,AS];
                    [~,ia,ic]=unique(tempPopulation.objs,'rows');
%                     childIndex=ic(1:size(Offspring,2));
                    tempPopulation=tempPopulation(ia);
%                     childIndex=Problem.N+1:length(tempPopulation);
                    tempZMax=zMax;
                    tempZMin=zMin;
                    index=find(tempZMax-tempZMin==0);
                    if(length(index)~=0)
                         tempZMax(index)=1;
                         tempZMin(index)=0;
                    end
                     range= tempZMax-tempZMin;
%                      fuzzy=5*10.^-5*(min(range,1));
                     epsilon=1*10.^-5*range;
                     fuzzy=epsilon;
%                      fuzzy=epsilon*(min(range,1));
%                      [FrontNo,MaxFNo] = NDSortTheta(tempPopulation.objs,min(Problem.N,length(tempPopulation)),fuzzy);
                    
                     [FrontNo,MaxFNo] = NDSort(ceil((tempPopulation.objs-zMin)./(fuzzy)),min(Problem.N,length(tempPopulation)));
%                      [FrontNo,MaxFNo] = NDSort((tempPopulation.objs),min(Problem.N,length(tempPopulation)));
%                      AddToAS=find(FrontNo(childIndex)==1);
%                      [FrontNo,MaxFNo] = NDSort(tempPopulation.objs,min(Problem.N,length(tempPopulation)));
%                      AddToAS=offspring(FrontNo(childIndex)==1);
                      

                      FristIndex=find(FrontNo==1);
                      normOfFirst=sum((tempPopulation(FristIndex).objs-zMin).^(2),2).^(1/2);
                      deleteIndex=normOfFirst>10*mean(normOfFirst);
                      FrontNo(FristIndex(deleteIndex))=2;
                      FristIndex(deleteIndex)=[];
%                       outlierIndex = determinedOutlier(tempPopulation(childIndex).objs-zMin,GNGnet);
%                       FrontNo(childIndex(outlierIndex))=max(FrontNo(childIndex(outlierIndex)),2);
                      [outlierIndex,minOfPoints,rate] = determinedOutlier(tempPopulation(FristIndex).objs,refPop,zMin,minOfPoints, alphaAngle,convexDgree,concaveDgree,rate);
%                       outlierIndex = determinedOutlier(tempPopulation.objs-zMin,AS.objs);  
                      FrontNo(FristIndex(outlierIndex))=max(FrontNo(FristIndex(outlierIndex)),2);

                      if(sum(FrontNo==1)>0)
                          AS=tempPopulation(FrontNo==1);
%                           tempFno=NDSort(AS.objs,1); % 没啥用。
%                           AS=AS(tempFno==1); % 没啥用。
                      else
                          
                          AS=tempPopulation(FrontNo==min(FrontNo));
%                           tempFno=NDSort(AS.objs,1); % 没啥用。
%                           AS=AS(tempFno==1); % 没啥用。
                      end
%                       if(AS==)
                      objs=AS.objs;
                      [zMax,~] = max(objs,[],1);
                      tempZMax=zMax;
                      tempZMin=zMin;
%                       try
                        index=find(tempZMax-tempZMin<10.^-20);
%                       catch
                        temp=0;
%                       end
                        if(length(index)~=0)
                            tempZMax(index)=1;
                            tempZMin(index)=0;
                        end
                      
                        
                      temp=(objs-tempZMin)./(tempZMax-tempZMin);
                      directionOfTemp=temp./sum(temp,2); % 以后不同方向都统一使用在超平面上的投影来进行计算。
%                       tempMax=max(temp,[],2);
                      indexOfM=1;
                      maxOfTemp=max(temp,[],2);
                      minOfTemp=min(temp,[],2);
%                       minOfTemp=zeros(1,Problem.M);
                      borderOfTemp=(maxOfTemp-minOfTemp);
                      cosineToLine=pdist2(directionOfTemp,ones(1,M)/M,'cosine');
%                       cosineToLine=pdist2(temp,ones(1,M),'cosine');
                      [~,maxIndex]=max(cosineToLine.*(temp(:,indexOfM)-minOfTemp));
                      cosineToLine(maxIndex)=[];
                      notInMaxIndex=setdiff(1:size(temp,1),maxIndex);
                      
%                       cosineToLine=cosineToLine./(temp(notInMaxIndex,indexOfM).^2);
                      indexOfM=indexOfM+1;
%                       cosineToLine=cosineToLine.*(temp(notInMaxIndex,indexOfM));
                      minCosDis=min(pdist2(directionOfTemp(notInMaxIndex,:),directionOfTemp(maxIndex(end),:),'cosine'),cosineToLine);
%                       minCosDis=min(pdist2(temp(notInMaxIndex,:),temp(maxIndex(end),:),'cosine'),cosineToLine);
                %       minCosDis=pdist2(temp(notInMaxIndex,:),temp(maxIndex(end),:),'cosine')+cosineToLine;

                      while length(maxIndex)<M-1 && length(notInMaxIndex)>0
%                           minCosDis= minCosDis.*(temp(notInMaxIndex,indexOfM).^2);
                          minCosDis=min(minCosDis,pdist2( directionOfTemp(notInMaxIndex,:), directionOfTemp(maxIndex(end),:),'cosine'));                       
                          [~,appendIndex]=max(minCosDis.*((temp(notInMaxIndex,indexOfM)- minOfTemp(notInMaxIndex))));
%                           [~,appendIndex]=max(minCosDis.*(tempMax(notInMaxIndex)));
                          maxIndex=[maxIndex,notInMaxIndex(appendIndex)];
                          minCosDis(appendIndex)=[];
                          notInMaxIndex(appendIndex)=[];
                          
%                           minCosDis=minCosDis./(temp(notInMaxIndex,indexOfM).^2);
                          indexOfM=indexOfM+1;
                      end
%                        %% 最后则侧重于夹角。
                       if length(maxIndex)<M && length(notInMaxIndex)>0
%                           minCosDis= minCosDis.*(temp(notInMaxIndex,indexOfM).^2);
                          minCosDis=min(minCosDis,pdist2( directionOfTemp(notInMaxIndex,:), directionOfTemp(maxIndex(end),:),'cosine'));                       
%                           [~,appendIndex]=max(minCosDis.*(temp(notInMaxIndex,indexOfM)-minOfTemp(notInMaxIndex)));
                          [~,appendIndex]=max(minCosDis.*(borderOfTemp(notInMaxIndex)));
                          maxIndex=[maxIndex,notInMaxIndex(appendIndex)];
                          minCosDis(appendIndex)=[];
                          notInMaxIndex(appendIndex)=[];
                          
%                           minCosDis=minCosDis./(temp(notInMaxIndex,indexOfM).^2);
%                           indexOfM=indexOfM+1;
                      end
%                      [~,maxIndex]=max(temp,[],1); % 传统的方法。

                     ExtremeVectors=temp(maxIndex,:);
%                    
                     
                     %% 方案3
                     q=0.03;
                     ExtremeVectors=ExtremeVectors-q*(mean(ExtremeVectors,1)-ExtremeVectors);
                     ExtremeVectors(ExtremeVectors>1)=1;
                     ExtremeVectors(ExtremeVectors<0)=0;
                     EM=ExtremeVectors;
                     EM=EM./(sum(EM,2));
                     %% 方案3


                     ChooseIndex = ArchiveUpdate_MD(AS.objs,ArchiveSize1,GNGnet,zMin,zMax,EM);
                     AS=AS(ChooseIndex);
                     
                     nAS = length(AS);
                     GNGnet.maxNode = min(ceil(1*Problem.N));
                     GNGnet.maxHP = 2*nAS; % paramter reset
                     GNGnet = GNGUpdateByDirection(AS.objs,GNGnet,Problem.M+1,zMin,ExtremeVectors);
%                      
                     %% Environmental Selection
                     [Population,FrontNo,GNGnet,zMax,crd] = ESelection_EM(tempPopulation,FrontNo,MaxFNo,Problem.N,Problem.M,Ru,GNGnet,zMin,zMax,0.75*(Problem.FE/(1*Problem.maxFE))+0.25,EM);
%                      GNGnet = GNGUpdateByDirection(Population.objs,GNGnet,Problem.M+1,zMin,ExtremeVectors);
                    
                else
%                     evaTimes=0;
%                     while(evaTimes<=Problem.N)
                        %% 原始方案。
%                         MatingPool = TournamentSelection(2,2,FrontNo,crd);
%             %           对整个优化使用小生境的方法来产生子代。
%                         try
%                             Offspring  = OperatorGAhalf(Problem,Population(MatingPool));
%                         catch
%                             temp=0;
%                         end
%                             evaTimes= evaTimes+1;
%                         zMin= min([zMin;Offspring.objs],[],1);
                        %% 此后使用基于LMPFE中拟合前沿的方式来依据投影点衡量个体的多样性以及收敛性。
                        % 估计前沿，首先利用GNG网络来判断前沿的分割区域，对不同的分割区域使用对应关联的种群个体来进行拟合。
                        if(ceil(Problem.FE/Problem.N) == MaxGen - NoG+1 && trainedFlag==0)
                            PopObjs=Population.objs;
                            FirstIndex=find(FrontNo==1);
                            objs=PopObjs(FirstIndex,:);
                            [zMax,maxIndex]=max(objs,[],1);
%                             ChooseIndex = ArchiveUpdate_MD(AS.objs,ArchiveSize1,GNGnet,zMin,zMax,EM);
%                             AS=AS(ChooseIndex);
%                             AS = ArchiveUpdate_MD([AS.objs;Population(FrontNo==1).objs],sum(FrontNo==1),ArchiveSize1,ArchiveSize1,GNGnet,zMin,zMax);
                             tempAS=AS;     
                            [Fronts,AS]=ApproximateFrontier(GNGnet,Population(FrontNo==1).objs,zMax,zMin);
%                             tempPopulation=[Population,AS];
% %                              tempPopulation=[Population,Offspring,AS];
%                             [~,ia,~]=unique(tempPopulation.objs,'rows');
%                             tempPopulation=tempPopulation(ia);
                            
%                             AS = ArchiveUpdate_MD([AS;Population(FrontNo==1).objs],sum(FrontNo==1),ArchiveSize1,ArchiveSize1,GNGnet,zMin,zMax);
%                             [tempFrontNo,~] = NDSortTheta(tempPopulation.objs,min(Problem.N,length(tempPopulation)),fuzzy);
%                             [Fronts,AS]=ApproximateFrontier(GNGnet,tempPopulation(tempFrontNo==1).objs,zMax,zMin);
%                           
%                             AS=tempPopulation(tempFrontNo==1);
%                             ChooseIndex = ArchiveUpdate_MD(AS.objs,ArchiveSize1,GNGnet,zMin,zMax,EM);
%                             AS=AS(ChooseIndex);
%                              tempAS=AS;       
%                             [Fronts,AS]=ApproximateFrontier(GNGnet,tempAS.objs,zMax,zMin);
                            %% 计算初始时的收敛度以及多样性的衡量指标。
                            [app,pp]=ProjectPoints(Fronts,Population.objs,AS,zMax,zMin);
                            [div,cov]=cpDivAndCov(Population.objs,app,zMin);
                            preDiv=mean(div);
                            preCov=mean(cov);
                            preDivs=[preDiv];
                            preCovs=[preCov];
                            theta2=0.7;
                            % theta2=abs(preCov/preDiv);
                            theta2s=[theta2];
                            trainedFlag=1;
                            
                            F1 = Population(FrontNo==1);
                            objs=F1.objs;
%                             [zMax,~]=max(objs,[],1);
                            [zMax,~] = max(AS,[],1);
                            info.zMax=zMax;
                            info.zMin=zMin;
                            info.TrueZMin=zMin;
                            
                            info.pp=pp;
                            PopObjs=(Population.objs-zMin)./(zMax-zMin);
                            info.cov=(app-sum(PopObjs,2))./(app);
                            info.dis=squareform(pdist(pp));
                            [sortDis,order]=sort(info.dis,2,'ascend');
                            info.sortDis=sortDis;
                            info.order=order;
                            info.Fno=FrontNo;
                            temp=0;
                            crd=zeros(1,Problem.N);
%                             [Population,FrontNo,info,crd,minOfPoints,rate]=ESelection_OneByOne1(tempAS,Problem,crd,Population,GNGnet,Fronts,info,minOfPoints, alphaAngle,convexDgree,concaveDgree,rate);
                            [Population,FrontNo,info,crd,minOfPoints,rate]=ESelection_OneByOneNew(tempAS,Problem,crd,Population,GNGnet,Fronts,info,minOfPoints, alphaAngle,convexDgree,concaveDgree,rate);
                        end
                      
                        %% 得到拟合前沿后，利用投影点进行筛选(但问题就在于估计不准，能与真实前沿的收敛程度达到0.1.....)，所以在由于还是使用基于支配等级优先。
                        % [Population,FrontNo,crd,EM]=ESelection_AF([Population,Offspring],Problem.N,Ruq,GNGnet,Fronts,NetLabel,theta2,zMin,EM);
                        % 版本5: 采用(miu+1)进化策略来选择。
                        zMin=info.zMin;
                        if(Problem.FE>=27073)
                            temp=0;
                        end
%                         [Population,FrontNo,info,crd]=ESelection_OneByOne(Problem,crd,Population,GNGnet,Fronts,info,epsilon);
                        MatingPool = TournamentSelection(2,Problem.N,info.Fno,crd);
                        Offsprings  = OperatorGA(Problem,Population(MatingPool)); % 直接一次性先产生N个个体。
                        [Population,FrontNo,info,crd,minOfPoints,rate]=ESelection_OneByOneNew(Offsprings,Problem,crd,Population,GNGnet,Fronts,info,minOfPoints, alphaAngle,convexDgree,concaveDgree,rate);
                        % 版本5: 采用(miu+1)进化策略来选择。

                        [app,~]=ProjectPoints(Fronts,Population.objs,AS,zMax,zMin);
                        [div,cov]=cpDivAndCov(Population.objs,app,zMin);


                        newDiv=mean(div);
                        newCov=mean(cov);
                        [theta2,preCov,preDiv]=UpdateTheta(preCov,preDiv,newCov,newDiv);
                        preDivs=[preDivs,preDiv];
                        preCovs=[preCovs,preCov];
                        theta2s=[theta2s;theta2];
%                     end
                end
            end
            temp=0;

            elapsed_time = toc;
%             fprintf('算法运行时间: %.4f 秒\n', elapsed_time);
            folder = fullfile('.\Data',class(Algorithm),class(Problem),num2str(Problem.M));
            [~,~]  = mkdir(folder);
            % folder=folder(1:end-1);
            file   = fullfile(folder,sprintf('elapsed_time_%s_%s_M%d_D%d_N%d_',class(Algorithm),class(Problem),Problem.M,Problem.D,Problem.N));
            runNo  = 1;
            while exist([file,num2str(runNo),'.mat'],'file') == 2
                runNo = runNo + 1;
            end
            save([file,num2str(runNo),'.mat'],'elapsed_time');            


            
            folder = fullfile('.\Data',class(Algorithm),class(Problem),num2str(Problem.M));
            [~,~]  = mkdir(folder);
            % folder=folder(1:end-1);
            file   = fullfile(folder,sprintf('GNGnetAndFronts_%s_%s_M%d_D%d_N%d_',class(Algorithm),class(Problem),Problem.M,Problem.D,Problem.N));
            runNo  = 1;
            while exist([file,num2str(runNo),'.mat'],'file') == 2
                runNo = runNo + 1;
            end
            nodes=GNGnet.NodeS;
            zMax=info.zMax;
            zMin=info.zMin;
            save([file,num2str(runNo),'.mat'],'GNGnet','tempAS','Fronts','preDivs','preCovs','theta2s','zMax','zMin');
            temp=0;
            
        end
    end
    
end