function [Centroids,AllPoints,RealClusterIDs] = mnl_GenerateDemoClusters_2d
Centroids=[1 2;1 8;2 2;3 8;4 1;4 5;4 6;7 1;8 4;9 7];
szCent=size(Centroids);
NumPointsPerClust=25;
AllPoints=nan(szCent(1)*NumPointsPerClust,szCent(2));
RealClusterIDs=nan(szCent(1)*NumPointsPerClust,1);
NoiseSpread=0.333;
for i=1:szCent
    tPoints=(randn(NumPointsPerClust,2).*NoiseSpread)+Centroids(i,:);
    %Squash points into the graph by making 0 and 10 the minimum and maximum respectively
    idx=tPoints<0;
    tPoints(idx)=0;
    idx=tPoints>10;
    tPoints(idx)=10;
    %Give a true cluster Id
    cID=ones(NumPointsPerClust,1)*i;
    %Allocate to group
    St=((i-1)*NumPointsPerClust)+1;
    Ed=St+NumPointsPerClust-1;
    AllPoints(St:Ed,:)=tPoints;
    RealClusterIDs(St:Ed,1)=cID;
    clear tPoints    
end
%% Plots with point colour
%Raw Values
nP=size(AllPoints,1);
Cmap=colormap(lines(szCent(1)));
figure('Name','Clusters by real grouping')
gscatter(AllPoints(:,1),AllPoints(:,2),RealClusterIDs,Cmap,'.',10)
grid on 
axis equal
xlim([0 10])
ylim([0 10])
xticks(0:10)
yticks(0:10)
legend('Location','northeastoutside')
end