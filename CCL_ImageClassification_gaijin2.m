% 单细胞和细胞簇分割开

function [Results] = CCL_ImageClassification_gaijin2(Mask, opts)
tic;

global imageName;

% 设置保存路径
path1 = 'C:\Users\lishaui\Desktop\CSC代码截图\path1';
path2 = 'C:\Users\lishaui\Desktop\CSC代码截图\path2';
path3 = 'C:\Users\lishaui\Desktop\CSC代码截图\path3';
path4 = 'C:\Users\lishaui\Desktop\CSC代码截图\singleSeed';  % 新增保存单一种子点区域的路径
path5 = 'C:\Users\lishaui\Desktop\CSC代码截图\multipleSeeds';  % 新增保存多种子点区域的路径
resultsFile = 'C:\Users\lishaui\Desktop\CSC代码截图\Results.txt'; % 结果保存的文本文件路径

threshold_low = 8;
threshold_high = 30;

CellsImage = zeros(size(Mask));

U = bwulterode(Mask);

[X, Y] = find(U > 0);

if (isfield(opts, 'max_num_seeds') == 0)
    opts.max_num_seeds = 1000;
end
if (isfield(opts, 'clp_th') == 0)
    opts.clp_th = 5;
end

if length(X) > opts.max_num_seeds
    Results.Labels = CellsImage;
    Results.nCells = 0;
    fprintf('Too many Seeds found: %d', length(X));
    return;
end

P = [X, Y];
[Seeds, Clusters, nclusters] = CSC_ClusterizePoints(P, opts.clp_th);

Q = zeros(size(Mask));
for i = 1:size(Seeds, 1)
    Q(Seeds(i, 1), Seeds(i, 2)) = i;
end

crem = bwmorph(Mask, 'remove');
CellSize = [];  % 存储细胞周长
CellArea = [];  % 存储细胞面积
t = [];
L = bwlabel(Mask > 0);
n = max(L(:));  % 连通区域数量

AreaList = [];  % 临时存储种子点数为1的区域的面积

% 新增区域分离部分
singleSeedCellsImage = zeros(size(Mask));  % 存储种子点数为1的区域
multipleSeedCellsImage = zeros(size(Mask));  % 存储种子点数大于1的区域

% 计算所有种子点数为1的区域的面积和周长
for k = 1:n
    labs = Q .* (L == k);
    nlabs = sum(labs(:) > 0);

    if (nlabs == 1) % 种子点数为1的区域
        clab = unique(labs(labs > 0));

        area = sum(L(:) == k); % 当前区域的面积
        perimeter = sum(sum(crem .* (L == k))); % 当前区域的周长

        % 将面积和周长分别添加到列表中
        AreaList = [AreaList; area];
        CellSize = [CellSize; perimeter];
        CellArea = [CellArea; area];
        t = [t, clab];
        CellsImage = CellsImage + clab .* (L == k);
        
        % 将该区域加入到单一种子点图像
        singleSeedCellsImage = singleSeedCellsImage + clab .* (L == k);
    elseif (nlabs > 1)  % 种子点数大于1的区域
    clab = unique(labs(labs > 0));

    area = sum(L(:) == k); % 当前区域的面积
    perimeter = sum(sum(crem .* (L == k))); % 当前区域的周长

    % 修正：对 `clab` 的值逐一更新 `multipleSeedCellsImage`
    for c = 1:length(clab)
        multipleSeedCellsImage = multipleSeedCellsImage + clab(c) .* (L == k);
    end
    end
end


% 对面积进行排序，并去掉最大10%和最小10%
[sortedAreas, sortIdx] = sort(AreaList);  % 对区域面积进行排序
N = length(sortedAreas);

if N > 10 % 确保区域数量足够多，避免过度筛选
    lowerIndex = ceil(0.1 * N); % 最小10%的索引
    upperIndex = floor(0.9 * N); % 最大10%的索引
    
    % 只保留中间80%的区域
    filteredIdx = sortIdx(lowerIndex:upperIndex); 
    filteredAreas = AreaList(filteredIdx);
    filteredSizes = CellSize(filteredIdx);
else
    filteredAreas = AreaList;
    filteredSizes = CellSize;
end

% 基于筛选后的数据重新计算统计值
AvgCellArea = mean(filteredAreas);      % 平均面积
StdCellArea = std(filteredAreas);       % 面积标准差
AvgCellSize = mean(filteredSizes);      % 平均周长
StdCellSize = std(filteredSizes);       % 周长标准差

% 打印统计结果
fprintf('AvgCellSize: %f\n', AvgCellSize);
fprintf('StdCellSize: %f\n', StdCellSize);
fprintf('AvgCellArea: %f\n', AvgCellArea);
fprintf('StdCellArea: %f\n', StdCellArea);

% 基于周长和面积计算细胞半径
Cell_r = (AvgCellSize + StdCellSize / 2) / (2 * pi);
Cellr_area = sqrt(AvgCellArea / pi);
fprintf('Cell_r: %f\n', Cell_r);
fprintf('Cellr_area: %f\n', Cellr_area);

% 将结果写入文件
fileID = fopen(resultsFile, 'a');
fprintf(fileID, '%s, AvgCellSize: %f, StdCellSize: %f, AvgCellArea: %f, StdCellArea: %f, Cellr_area: %f, Cell_r: %f\n\n', ...
    imageName, AvgCellSize, StdCellSize, AvgCellArea, StdCellArea, Cellr_area, Cell_r);
fclose(fileID);

% 保存种子点数为1的区域图像
imwrite(singleSeedCellsImage, fullfile(path4, sprintf('%s_singleSeed.png', imageName)));
fprintf('Single seed image saved to path4.\n');

% 保存种子点数大于1的区域图像
imwrite(multipleSeedCellsImage, fullfile(path5, sprintf('%s_multipleSeeds.png', imageName)));
fprintf('Multiple seeds image saved to path5.\n');

% 保存图像到相应路径
if StdCellSize > threshold_high
    imwrite(Mask, fullfile(path1, sprintf('%s_%f.png', imageName, StdCellSize)));
    fprintf('Image saved to path1.\n');
elseif StdCellSize < threshold_low
    imwrite(Mask, fullfile(path3, sprintf('%s_%f.png', imageName, StdCellSize)));
    fprintf('Image saved to path3.\n');
else
    imwrite(Mask, fullfile(path2, sprintf('%s_%f.png', imageName, StdCellSize)));
    fprintf('Image saved to path2.\n');
end

% 输出结果结构体
Results.AvgCellSize = AvgCellSize;
Results.StdCellSize = StdCellSize;
Results.AvgCellArea = AvgCellArea;
Results.StdCellArea = StdCellArea;
Results.cell_r = Cell_r;
Results.Labels = CellsImage;
Results.Time = toc;

end

%% %%%%%%%%%%%%%%%%%% Side Functions辅助功能 %%%%%%%%%%%%%%%%%%%%%%%%%%% %%

function [Seeds, Clusters, nclusters] = CSC_ClusterizePoints(P, th)
    % 计算种子点之间的距离矩阵 M。
    n = max(size(P));

    if(n > 1)
        M = zeros(n, n);
        for i=1:n
            for j=1:n
                M(i,j)=abs(P(i,1)-P(j,1))+abs(P(i,2)-P(j,2));
            end
        end
    end

    % 对点进行聚类，如果点之间的距离小于 th，则将它们归为同一类，否则创建新的聚类。
    Clusters = zeros(1,n);
    dist = zeros(1,n);
    nclusters = 0;
    for i=1:n
        dist(:)=Inf;
        for j=1:nclusters
            k=1;
            while(Clusters(j,k)>0)
                dist(j) = min(dist(j), M(i,Clusters(j,k)));
                k = k+1;
            end
        end
        dminpos = find(dist==min(dist),1);
        dmin = min(dist);
        if(dmin <= th)
            pos = find(Clusters(dminpos,:)==0,1);
            Clusters(dminpos, pos) = i;
        else
            nclusters = nclusters + 1;
            Clusters(nclusters, 1) = i;
        end
    end

    % 计算每个聚类的质心作为种子点，返回这些种子点的坐标。
    Seeds = zeros(nclusters,2);
    for i=1:nclusters
        pos = find(Clusters(i,:)==0,1)-1;
        Seeds(i,1) = round(mean(P(Clusters(i,1:pos),1)));
        Seeds(i,2) = round(mean(P(Clusters(i,1:pos),2)));
    end
end
