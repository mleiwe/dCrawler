function [ClusterIDs,Centroids]=mnl_Basic_dCrawler(InputMatrix,EuThresh,FigYN)
% Function to create clusters based on the Euclidean Distances, this is the
% basic form with just a matrix provided. Now with improved speed on the
% Adjust clusters section. Then after the merge clusters the points outside
% the cluster are now re-allocated
% 
% Inputs
%  InputMatrix - the list of points, each row = point, column = values
%  EuThresh - the euclidean threshold (d) you want to use
%  FigYN - do you want figures or not ('y'/'n')
% 
% Outputs
%  ClusterIDs - The id of each trace in a single matrix, each row is the points cluster
%  Centroids - The centroid position of each cluster
% 
% This version has an improved video quality output
% Pre-allocation
i=1;
ClusterNum=1;
nCol=size(InputMatrix,2);
nRow=size(InputMatrix,1);
IntermediateMatrix=InputMatrix;
Points=[];
ClusterIDs=zeros(nRow,1);
% Asign the first values
NewPoints=InputMatrix(i,:);
IntermediateMatrix(i,:)=nan(1,nCol);
ClusterIDs(i)=ClusterNum;
Points=[Points;NewPoints];
Centroids(i,:)=Points;
TracesPerCluster=1;
TracesOfCluster(TracesPerCluster)=i;
% Check if the user wants a figure
if strcmp(FigYN,'y')==1
    prompt=sprintf('%s','Do you want to view intermediate steps? NB This will slow down the code alot (y/n)');
    FigInt=input(prompt,'s');
end
if strcmp(FigYN,'y')==1
    M=VideoWriter('DemoMovie.avi');
    M.FrameRate=10;   
    f=1;
    open(M)
    mkdir('Movie_JPGs')
    if nCol>2 %if nCol>2, then tSNE it
        Y=tsne(InputMatrix);
    else
        Y=InputMatrix;
    end
    minX=floor(min(Y(:,1)));
    maxX=ceil(max(Y(:,1)));
    minY=floor(min(Y(:,2)));
    maxY=ceil(max(Y(:,2)));
    minV=min([minX minY]);
    maxV=max([maxX maxY]);
    figure('Name','Clusters','Units','normalized','Position',[0.1 0.1 0.33 0.5])
    scatter(Y(:,1),Y(:,2),'.k')
    hold on
    grid on
    axis equal
    xlim([minV maxV])
    ylim([minV maxV])
    xticks(minV:maxV)
    yticks(minV:maxV)
    legend('Location','northeastoutside')
    %Now save for movie
    set(0,'defaultlinelinesmoothing');
    fname=['Movie_JPGs\MovieFrame' num2str(f)];
    print(fname,'-djpeg','-r200');
    I=imread([fname '.jpg']);
    frame=im2frame(I);
    writeVideo(M,frame)
end    
% Step 1 - The first Crawl through
while i<=nRow
    Centroid=mean(Points,1,'omitnan');
    Centroids(ClusterNum,:)=Centroid;
    if strcmp(FigYN,'y')==1
        if strcmp(FigInt,'y')==1
            %Create the colormap
            tCmap=colormap(lines(ClusterNum));
            Cmap(1,:)=[0 0 0]; %Make the first (unassigned) group black
            Cmap(2:ClusterNum+1,:)=tCmap;
            %Plot points
            colormap(Cmap)
            [ii,~]=find(~ClusterIDs);
            if isempty(ii)==1
                gscatter(Y(:,1),Y(:,2),ClusterIDs,tCmap,'.',15)
            else
                gscatter(Y(:,1),Y(:,2),ClusterIDs,Cmap,'.',15)
            end
            tn=sprintf('%s%d','Initial Crawl...Cluster ',ClusterNum);
            title(tn)
            if nCol==2
                if exist('hCent','var')==1
                    delete(hCent)
                    delete(hCirc)
                end
                %Plot Centroids
                legend off
                hCent=gscatter(Centroids(:,1),Centroids(:,2),1:ClusterNum,tCmap,'x',20,'off');
                %Plot Radii
                for j=1:ClusterNum
                    hCirc(j)=viscircles(Centroids(j,:),EuThresh,'LineStyle','--','Color',tCmap(j,:));
                end 
            end
            xlim([minV maxV]);
            ylim([minV maxV]);
            for name=1:ClusterNum
                ClustNames{name}=num2str(name);
            end
            legend(hCent,ClustNames,'Location','northeastoutside')
            %Now save for movie
            f=f+1;
            set(0,'defaultlinelinesmoothing');
            fname=['Movie_JPGs\MovieFrame' num2str(f)];
            print(fname,'-djpeg','-r200');
            I=imread([fname '.jpg']);
            frame=im2frame(I);
            writeVideo(M,frame);
        end
    end
    
    % Check if the points are still close enough to the centroid
    [tDistMatrix]=mnl_GroupEuD_v2(Centroid,Points); %This is faster
    [Pos]=find(tDistMatrix>EuThresh);
    if isempty(Pos)==0
        NumberOutside=size(Pos,2);
        %Return Values back for assessment
        for j=1:NumberOutside
            %Get the Row Id
            RowID=TracesOfCluster(Pos(j));
            %Return the colour to the comparing matrix
            IntermediateMatrix(RowID,:)=InputMatrix(RowID,:);
            %Return the cluster ids
            ClusterIDs(RowID)=0;
            %Now Mark the points for deletion
            Points(Pos(j),:)=nan(1,nCol); %initially mark them NaN
            TracesOfCluster(Pos(j))=NaN;
        end
        %Remove NaNs and reset the counters
        szPoints=size(Points,1);
        n=1;
        nC=1;
        for j=1:szPoints
            if ~isnan(Points(j,:))
                tPoints(n,:)=Points(j,:);
                n=n+1;
            end
            if ~isnan(TracesOfCluster(j))
                tTracesOfCluster(nC)=TracesOfCluster(j);
                nC=nC+1;
            end
        end
        Points=tPoints;
        TracesOfCluster=tTracesOfCluster;
        TracesPerCluster=TracesPerCluster-NumberOutside;
        %Now Recalculate the centroid
        Centroid=mean(Points,1,'omitnan');
        clear Pos tPoints tTracesOfCluster
    end    
    %Calculate the distance to the Center
    [DistMatrix]=mnl_GroupEuD_v2(Centroid,IntermediateMatrix); %Faster
    %Find the minumum distance
    [val,loc]=min(DistMatrix);
    %Is it below the Euclidean Distance
    if val<=EuThresh
        %Update the location
        ClusterIDs(loc)=ClusterNum;
        NewPoints=IntermediateMatrix(loc,:);
        Points=[Points;NewPoints];
        IntermediateMatrix(loc,:)=nan(1,nCol);       
        TracesPerCluster=TracesPerCluster+1;
        TracesOfCluster(TracesPerCluster)=loc;
    else
        %Save as a cluster
        Centroids(ClusterNum,:)=Centroid;
        clear TracesOfCluster NewPoints
        %Make A new Cluster
        ClusterNum=ClusterNum+1;
        TracesPerCluster=1;
        %Move to the next value
        i=1;
        while 1            
            if i<=nRow
                if ~isnan(IntermediateMatrix(i,:)) %keep cycling through until we find a value that hasn't been allocated to a cluster
                    [a,~]=find(ClusterIDs~=0);
                    numDefined=size(a,1);
                    mnl_InsertProgressTrackerInLoops(numDefined,nRow)                   
                    break
                else
                    i=i+1;
                end
            else
                break
            end
        end
        %Redefine if i<=nRow
        if i<=nRow
            Points=[];
            NewPoints=IntermediateMatrix(i,:);
            IntermediateMatrix(i,:)=nan(1,nCol);
            ClusterIDs(i)=ClusterNum;
            Points=[Points;NewPoints];
            TracesPerCluster=1;
            TracesOfCluster(1)=i;
        else
            break
        end
    end    
end
% Ammend it to make sure the points are associated to their nearest cluster
ClusterNum=ClusterNum-1;
i=1;
disp('Adjusting Clusters...')
LoopAdjust=0;
nLoop=0;
while i<=nRow %NB Cannot just loop it need to keep checking at each stage until all are associated with their nearest centroid
    nAdjust=0;
    Pos1=InputMatrix(i,:);
    Clus=ClusterIDs(i);
    [DistMatrix]=mnl_GroupEuD_v2(Pos1,Centroids);
    [~,loc]=min(DistMatrix);
    if Clus~=loc
        %Re-assign ClusterID
        ClusterIDs(i)=loc;
        %Re-calculate the new centroid
        [RowList]=find(ClusterIDs==loc);
        tPoints=InputMatrix(RowList,:);
        Centroids(loc,:)=mean(tPoints,1);        
        %Re-calculate the old centroid
        [RowList]=find(ClusterIDs==Clus);
        tPoints=InputMatrix(RowList,:);
        Centroids(Clus,:)=mean(tPoints,1);
        %Update the counters
        nAdjust=nAdjust+1;
        LoopAdjust=LoopAdjust+1;
        if strcmp(FigYN,'y')==1
            if strcmp(FigInt,'y')==1
                tCmap=colormap(lines(ClusterNum));
                Cmap(1,:)=[0 0 0]; %Make the first (unassigned) group black
                Cmap(2:ClusterNum+1,:)=tCmap;
                tn=sprintf('%s','Adjusting Clusters...Allocate each point to their nearest centroid');
                title(tn)
                gscatter(Y(:,1),Y(:,2),ClusterIDs,tCmap,'.',15)
                hold on
                if nCol==2
                    if exist('hCent','var')==1
                        delete(hCent)
                        delete(hCirc)
                    end
                    %Plot Centroids
                    legend off
                    hCent=gscatter(Centroids(:,1),Centroids(:,2),1:ClusterNum,tCmap,'x',20,'off');
                    %Plot Radii
                    for j=1:ClusterNum
                        hCirc(j)=viscircles(Centroids(j,:),EuThresh,'LineStyle','--','Color',tCmap(j,:));
                    end
                end
                xlim([minV maxV]);
                ylim([minV maxV]);
                for name=1:ClusterNum
                    ClustNames{name}=num2str(name);
                end
                legend(hCent,ClustNames,'Location','northeastoutside')
                %Now save for movie
                f=f+1;
                set(0,'defaultlinelinesmoothing');
                fname=['Movie_JPGs\MovieFrame' num2str(f)];
                print(fname,'-djpeg','-r200');
                I=imread([fname '.jpg']);
                frame=im2frame(I);
                writeVideo(M,frame);
            end
        end
        if nAdjust>0
            if i+1<nRow
                i=i+1; %inserted so we loop back to the begining
            else
                i=1;
            end
        end
    elseif Clus==loc %Then move on to the next point
        if i+1<nRow
            i=i+1;
        else %if we have reached the end of the loop
            if LoopAdjust>0 %If we found points to adjust 
                nLoop=nLoop+1;
                str=sprintf('%s%d%s%d%s','Done ',nLoop,' Loops with ',LoopAdjust,' Adjustments');
                disp(str)
                LoopAdjust=0;
                i=1;
            else
                i=i+1;
            end
        end
    end
end
% Merge Clusters if they are too close
ClusterThresh=EuThresh;
if strcmp(FigYN,'y')==1
    PermCmap=colormap(lines(ClusterNum));
end
disp('Merging Nearby Clusters')
mxClusterNum=ClusterNum;
%Measure the distances between clusters
mnLoop=1;
while 1
    nMerged=0;
    i=1;
    while i<=mxClusterNum
        Pos1=Centroids(i,:);
        tCent=Centroids;
        tCent(i,:)=nan(1,nCol);
        [DistMatrix]=mnl_GroupEuD_v2(Pos1,tCent);
        [loc]=find(DistMatrix<ClusterThresh);
        if isempty(loc)~=1
            nMerged=nMerged+1;
            %Now merge the closest one
            [MinVal,MinLoc]=min(DistMatrix(loc));
            PotentialMerge=[i;loc(MinLoc)];
            NewClusterNumber=min(PotentialMerge);
            %Calculate the new centroid position
            tClustPoints=[];
            for j=1:2
                tidx=ClusterIDs==PotentialMerge(j);
                tClustPoints=[tClustPoints;InputMatrix(tidx,:)];
            end
            TempCentroidPos=mean(tClustPoints);
            clear tClustPoints
            %Now recode the clusters
            for j=1:2
                if NewClusterNumber==PotentialMerge(j)
                    %Is the point close enough to the centroid?
                    [RowList]=find(ClusterIDs==PotentialMerge(j));
                    tMatrix=InputMatrix(RowList,:);
                    [DistMatrix]=mnl_GroupEuD_v2(TempCentroidPos,tMatrix);
                    %Allocate those within with the new number
                    [locInc]=find(DistMatrix<=EuThresh);
                    [locExc]=find(DistMatrix>EuThresh);
                    if isempty(locInc)==0
                        ClusterIDs(RowList(locInc))=NewClusterNumber;
                    end
                    %Score those outside as zero
                    if isempty(locExc)==0
                        ClusterIDs(RowList(locExc))=0;
                    end
                    clear tMatrix
                elseif NewClusterNumber~=PotentialMerge(j)
                    [RowList]=find(ClusterIDs==PotentialMerge(j));
                    tMatrix=InputMatrix(RowList,:);
                    [DistMatrix]=mnl_GroupEuD_v2(TempCentroidPos,tMatrix);
                    %Allocate those within with the new number
                    [locInc]=find(DistMatrix<=EuThresh);
                    [locExc]=find(DistMatrix>EuThresh);
                    if isempty(locInc)==0
                        ClusterIDs(RowList(locInc))=NewClusterNumber;
                    end
                    %Score those outside as zero
                    if isempty(locExc)==0
                        ClusterIDs(RowList(locExc))=0;
                    end
                    clear tMatrix
                    %Score the discarded centroids as nan
                    Centroids(PotentialMerge(j),:)=nan(1,nCol);
                end
            end
            clear TempCentroidPos
            %Update the centroids
            %Re-calculate the new centroid position
            [RowList]=find(ClusterIDs==NewClusterNumber);
            tPoints=InputMatrix(RowList,:);
            Centroids(NewClusterNumber,:)=mean(tPoints,1);
            clear tPoints           
            %Update the figure
            if strcmp(FigYN,'y')==1
                if strcmp(FigInt,'y')==1
                    %Re-code the colourmap for points
                    g=1;
                    gC=1;
                    PointMap=[];
                    CentMap=[];
                    for k=1:mxClusterNum
                        if isnan(Centroids(k,:))~=1
                            PointMap(g,:)=PermCmap(k,:);
                            CentMap(gC,:)=PermCmap(k,:);
                            g=g+1;
                            gC=gC+1;
                        elseif isnan(Centroids(k,:))>=1
                            CentMap(gC,:)=[0 0 0];
                            gC=gC+1;
                        end
                    end
                    tCmap=PointMap;
                    if min(ClusterIDs)<1
                        Cmap(1,:)=[0 0 0]; %Make the first (unassigned) group black
                        Cmap(2:g,:)=tCmap;
                    else
                        Cmap=tCmap;
                    end
                    tn=sprintf('%s%d%s%d','Merging Clusters...',PotentialMerge(1),'and ',PotentialMerge(2));
                    title(tn)
                    gscatter(Y(:,1),Y(:,2),ClusterIDs,Cmap,'.',15)
                    hold on
                    if nCol==2
                        if exist('hCent','var')==1
                            delete(hCent)
                            delete(hCirc)
                        end
                        %Plot Centroids
                        legend off
                        hCent=gscatter(Centroids(:,1),Centroids(:,2),1:mxClusterNum,CentMap,'x',15,'off');
                        %Plot Radii
                        for j=1:mxClusterNum
                            hCirc(j)=viscircles(Centroids(j,:),EuThresh,'LineStyle','--','Color',CentMap(j,:));
                        end
                    end
                    xlim([minV maxV]);
                    ylim([minV maxV]);
                    for name=1:ClusterNum
                        ClustNames{name}=num2str(name);
                    end
                    legend(hCent,ClustNames,'Location','northeastoutside')
                    %Now save for movie
                    f=f+1;
                    set(0,'defaultlinelinesmoothing');
                    fname=['Movie_JPGs\MovieFrame' num2str(f)];
                    print(fname,'-djpeg','-r200');
                    I=imread([fname '.jpg']);
                    frame=im2frame(I);
                    writeVideo(M,frame);
                end
            end
        end  
        i=i+1;
    end
    %Are There any Points that are unlabelled?
    [RowList]=find(ClusterIDs==0);
    %If yes then we need to crawl through again
    if isempty(RowList)==0
        %Find the new cluster number
        ClustNumList=unique(ClusterIDs,'rows');
        BiggestNumber=max(ClusterIDs);
        c=1;
        for j=1:BiggestNumber
            [loc]=find(ClustNumList==j);
            if isempty(loc)==1
                MissingClustNum(c)=j;
                c=c+1;
            end
        end
        nMissingClusters=c-1;
        %Create a New Matrix for unassigned points
        TempPoints=InputMatrix(RowList,:);
        nPoints=size(TempPoints,1);
        %Make new clusters
        if nMissingClusters>0
            ClusterNum=MissingClustNum(1);
        else
            ClusterNum=max(ClustNumList)+1;
        end
        %Pre-allocate
        j=1;
        Points=[];
        NewPoints=TempPoints(j,:);
        IntermediateMatrix(RowList,:)=TempPoints;
        Points=[Points;NewPoints];
        IntermediateMatrix(RowList(j),:)=nan(1,nCol);
        ClusterIDs(RowList(j))=ClusterNum;
        TracesOfCluster(TracesPerCluster)=RowList(j);
        TracesPerCluster=1;
        %Now Explore the unassigned for crawling
        while j<=nRow
            Centroid=mean(Points,1,'omitnan');
            Centroids(ClusterNum,:)=Centroid;
            if strcmp(FigYN,'y')==1
                if strcmp(FigInt,'y')==1
                    %Update the colormap
                    tCmap=colormap(lines(mxClusterNum));
                    g=1;
                    gC=1;
                    PointMap=[];
                    CentMap=[];
                    for k=1:mxClusterNum
                        if isnan(Centroids(k,:))~=1
                            PointMap(g,:)=PermCmap(k,:);
                            CentMap(gC,:)=PermCmap(k,:);
                            g=g+1;
                            gC=gC+1;
                        elseif isnan(Centroids(k,:))>=1
                            CentMap(gC,:)=[0 0 0];
                            gC=gC+1;
                        end
                    end
                    tCmap=PointMap;
                    if min(ClusterIDs)<1 %If there are still unassigned points
                        Cmap(1,:)=[0 0 0]; %Make the first (unassigned) group black
                        Cmap(2:g,:)=tCmap;
                    else
                        Cmap=tCmap;
                    end
                    %Plot points
                    colormap(tCmap)
                    [ii,~]=find(~ClusterIDs);
                    if isempty(ii)==1
                        gscatter(Y(:,1),Y(:,2),ClusterIDs,tCmap,'.',15)
                    else
                        gscatter(Y(:,1),Y(:,2),ClusterIDs,Cmap,'.',15)
                    end
                    tn=sprintf('%s%d','Reassigning Points Outside of Cluster ',ClusterNum);
                    title(tn)
                    if nCol==2
                        if exist('hCent','var')==1
                            delete(hCent)
                            delete(hCirc)
                        end
                        %Plot Centroids
                        legend off
                        hCent=gscatter(Centroids(:,1),Centroids(:,2),1:mxClusterNum,CentMap,'x',20,'off');
                        %Plot Radii
                        for k=1:mxClusterNum
                            hCirc(k)=viscircles(Centroids(k,:),EuThresh,'LineStyle','--','Color',CentMap(k,:));
                        end
                    end
                    xlim([minV maxV]);
                    ylim([minV maxV]);
                    for name=1:ClusterNum
                        ClustNames{name}=num2str(name);
                    end
                    legend(hCent,ClustNames,'Location','northeastoutside')
                    %Now save for movie
                    f=f+1;
                    set(0,'defaultlinelinesmoothing');
                    fname=['Movie_JPGs\MovieFrame' num2str(f)];
                    print(fname,'-djpeg','-r200');
                    I=imread([fname '.jpg']);
                    frame=im2frame(I);
                    writeVideo(M,frame);
                end
            end
            [tDistMatrix]=mnl_GroupEuD_v2(Centroid,Points);
            [Pos]=find(tDistMatrix>EuThresh);
            if isempty(Pos)==0
                NumberOutside=size(Pos,2);
                %Return Values back for assessment
                for k=1:NumberOutside
                    %Get the Row Id
                    RowID=TracesOfCluster(Pos(k));
                    %Return the colour to the comparing matrix
                    IntermediateMatrix(RowID,:)=InputMatrix(RowID,:);
                    %Return the cluster ids
                    ClusterIDs(RowID)=0;
                    %Now Mark the points for deletion
                    Points(Pos(k),:)=nan(1,nCol); %initially mark them NaN
                    TracesOfCluster(Pos(k))=NaN;
                end
                %Remove NaNs and reset the counters
                szPoints=size(Points,1);
                n=1;
                nC=1;
                for k=1:szPoints
                    if isnan(Points(k,:))==0
                        tPoints(n,:)=Points(k,:);
                        n=n+1;
                    end
                    if isnan(TracesOfCluster(k))==0
                        tTracesOfCluster(nC)=TracesOfCluster(k);
                        nC=nC+1;
                    end
                end
                Points=tPoints;
                TracesOfCluster=tTracesOfCluster;
                TracesPerCluster=TracesPerCluster-NumberOutside;
                %Now Recalculate the centroid
                Centroid=nanmean(Points,1);
                clear Pos tPoints tTracesOfCluster
            end
            %Calculate the distance to the Center
            [DistMatrix]=mnl_GroupEuD_v2(Centroid,IntermediateMatrix);
            %Find the minumum distance
            [val,loc]=min(DistMatrix);
            %Is it below the Euclidean Distance
            if val<=EuThresh
                %Update the location
                ClusterIDs(loc)=ClusterNum;
                NewPoints=IntermediateMatrix(loc,:);
                Points=[Points;NewPoints];
                IntermediateMatrix(loc,:)=nan(1,nCol);
                TracesPerCluster=TracesPerCluster+1;
                TracesOfCluster(TracesPerCluster)=loc;
            else
                %Save as a cluster
                Centroids(ClusterNum,:)=Centroid;
                clear TracesOfCluster NewPoints
                %Make A new Cluster
                %Find the new cluster number
                ClustNumList=unique(ClusterIDs,'rows');
                BiggestNumber=max(ClusterIDs);
                c=1;
                MissingClustNum=[];
                for k=1:BiggestNumber
                    [loc]=find(ClustNumList==k);
                    if isempty(loc)==1
                        MissingClustNum(c)=k;
                        c=c+1;
                    end
                end
                nMissingClusters=c-1;
                if c>1
                    ClusterNum=min(MissingClustNum);
                else               
                    ClusterNum=BiggestNumber+1;
                end
                TracesPerCluster=1;
                %Move to the next value
                j=1;
                while 1
                    if j<=nRow
                        if isnan(IntermediateMatrix(j,:))==0 %keep cycling through until we find a value that hasn't been allocated to a cluster
                            [a,~]=find(ClusterIDs~=0);
                            numDefined=size(a,1);
                            mnl_InsertProgressTrackerInLoops(numDefined,nRow)
                            break
                        else
                            j=j+1;
                        end
                    else
                        break
                    end
                end
                %Redefine if j<=nRow
                if j<=nRow
                    Points=[];
                    NewPoints=IntermediateMatrix(j,:);
                    IntermediateMatrix(j,:)=nan(1,nCol);
                    ClusterIDs(j)=ClusterNum;
                    Points=[Points;NewPoints];
                    TracesPerCluster=1;
                    TracesOfCluster(1)=j;
                else
                    break
                end
            end
        end
        %Adjust the new clusters
        j=1;
        disp('Readjusting Clusters After Merging...')
        LoopAdjust=0;
        nLoop=0;
        while j<=nRow %NB Cannot just loop it need to keep checking at each stage until all are associated with their nearest centroid
            nAdjust=0;
            Pos1=InputMatrix(j,:);
            Clus=ClusterIDs(j);
            %[DistMatrix]=mnl_GroupEuD(Pos1,Centroids);
            [DistMatrix]=mnl_GroupEuD_v2(Pos1,Centroids);
            [~,loc]=min(DistMatrix);
            mxClusterNum=max(ClusterIDs);%redefine the max cluster number
            if Clus~=loc
                %Re-assign ClusterID
                ClusterIDs(j)=loc;
                %Re-calculate the new centroid
                [RowList]=find(ClusterIDs==loc);
                tPoints=InputMatrix(RowList,:);
                Centroids(loc,:)=mean(tPoints,1);
                %Re-calculate the old centroid
                [RowList]=find(ClusterIDs==Clus);
                tPoints=InputMatrix(RowList,:);
                Centroids(Clus,:)=mean(tPoints,1);
                %Update the counters
                nAdjust=nAdjust+1;
                LoopAdjust=LoopAdjust+1;
                if strcmp(FigYN,'y')==1
                    if strcmp(FigInt,'y')==1
                        %Update the colormap
                        tCmap=colormap(lines(mxClusterNum));
                        g=1;
                        gC=1;
                        PointMap=[];
                        CentMap=[];
                        for k=1:mxClusterNum
                            if isnan(Centroids(k,:))~=1
                                PointMap(g,:)=PermCmap(k,:);
                                CentMap(gC,:)=PermCmap(k,:);
                                g=g+1;
                                gC=gC+1;
                            elseif isnan(Centroids(k,:))>=1
                                CentMap(gC,:)=[0 0 0];
                                gC=gC+1;
                            end
                        end
                        tCmap=PointMap;
                        if min(ClusterIDs)<1 %If there are still unassigned points
                            Cmap(1,:)=[0 0 0]; %Make the first (unassigned) group black
                            Cmap(2:g,:)=tCmap;
                        else
                            Cmap=tCmap;
                        end
                        tn=sprintf('%s','Readjusting Clusters After Merging...Allocating each point to their nearest centroid');
                        title(tn)
                        gscatter(Y(:,1),Y(:,2),ClusterIDs,PointMap,'.',15)
                        hold on
                        if nCol==2
                            if exist('hCent','var')==1
                                delete(hCent)
                                delete(hCirc)
                            end
                            %Plot Centroids
                            legend off
                            hCent=gscatter(Centroids(:,1),Centroids(:,2),1:mxClusterNum,CentMap,'x',20,'off');
                            %Plot Radii
                            for j=1:mxClusterNum
                                hCirc(j)=viscircles(Centroids(j,:),EuThresh,'LineStyle','--','Color',CentMap(j,:));
                            end
                        end
                        xlim([minV maxV]);
                        ylim([minV maxV]);
                        for name=1:ClusterNum
                            ClustNames{name}=num2str(name);
                        end
                        legend(hCent,ClustNames,'Location','northeastoutside')
                        %Now save for movie
                        f=f+1;
                        set(0,'defaultlinelinesmoothing');
                        fname=['Movie_JPGs\MovieFrame' num2str(f)];
                        print(fname,'-djpeg','-r200');
                        I=imread([fname '.jpg']);
                        frame=im2frame(I);
                        writeVideo(M,frame);
                    end
                end
                if nAdjust>0
                    if j+1<nRow
                        j=j+1; %move on to the next point
                    else
                        j=1;%inserted so we loop back to the begining
                    end
                end
            elseif Clus==loc %It is correct so move on to the next point
                if j+1<nRow
                    j=j+1; %move on to the next point
                else %if we have reached the end of the loop
                    if LoopAdjust>0 %If we found points to adjust
                        nLoop=nLoop+1;
                        str=sprintf('%s%d%s%d%s','Done ',nLoop,' Loops with ',LoopAdjust,' Adjustments');
                        disp(str)
                        LoopAdjust=0;
                        j=1;
                    else
                        j=j+1;
                    end
                end
            end
        end
    end
    %Now check if another loop is needed
    if nMerged>=1
        fprintf('%s%d%s%d%s','Merge Loop ',mnLoop,' Merged ',nMerged,' Cluster Sets');
        mnLoop=mnLoop+1;
        mxClusterNum=max(ClusterIDs);
        i=1; %Now Start from the begining 
        %if it has looped more than 20 times reduce the merging threshold by 1%
        if mnLoop>=21
            ClusterThresh=ClusterThresh*0.95;
        end
    elseif nMerged==0
        break
    end
end
% Final Figure
%Clean up the cluster numbers
disp('Cleaning up the Cluster Numbers')
mxClusterNum=max(ClusterIDs);
AllNums=unique(ClusterIDs,'rows');
for i=1:mxClusterNum
    %Does it exist already?
    loc=find(AllNums==i);
    if isempty(loc)==1 && i<=max(AllNums)
        %Allocate the next bigest number
        idx=AllNums(:)>i;
        [~,Pos]=max(idx);
        OldClusterNum=AllNums(Pos);
        %Update Centroids
        Centroids(i,:)=Centroids(OldClusterNum,:);
        Centroids(OldClusterNum,:)=nan;
        %Update ClusterIDs
        [RowList]=find(ClusterIDs==OldClusterNum);
        ClusterIDs(RowList)=i;
        %Update AllNums
        AllNums=unique(ClusterIDs,'rows');
    end
end
[r,c]=find(isnan(Centroids)==0);
Centroids2(r,c)=Centroids(r,c);
Centroids=Centroids2; clear Centroids2
if strcmp(FigYN,'y')==1
    %Re-code the colourmap for points - This is the old code where the ClusterIDs and Centroids may have contained nans
    tCmap=colormap(lines(mxClusterNum));
    g=1;
    gC=1;
    PointMap=[];
    CentMap=[];
    szfCentroids=size(Centroids,1);
    for k=1:szfCentroids
        if ~isnan(Centroids(k,:))
            PointMap(g,:)=PermCmap(k,:);
            CentMap(gC,:)=PermCmap(k,:);
            g=g+1;
            gC=gC+1;
        elseif isnan(Centroids(k,:))
            CentMap(gC,:)=[0 0 0];
            gC=gC+1;
        end
    end
    tCmap=PointMap;
    mxClusterNum=max(ClusterIDs);
    PointMap=colormap(lines(mxClusterNum));                
    gscatter(Y(:,1),Y(:,2),ClusterIDs,PointMap,'.',15)
    title('Final Clusters')
    hold on
    if nCol==2
        if exist('hCent','var')==1
            delete(hCent)
            delete(hCirc)
        end
        %Plot Centroids
        hCent=gscatter(Centroids(:,1),Centroids(:,2),1:mxClusterNum,PointMap,'x',15);
        %Plot Radii
        for j=1:mxClusterNum
            hCirc(j)=viscircles(Centroids(j,:),EuThresh,'LineStyle','--','Color',PointMap(j,:));
        end
    end
    xlim([minV maxV]);
    ylim([minV maxV]);
    for name=1:ClusterNum
        ClustNames{name}=num2str(name);
    end
    legend(hCent,ClustNames,'Location','northeastoutside')
    %Now save for movie
    f=f+1;
    set(0,'defaultlinelinesmoothing');
    fname=['Movie_JPGs\MovieFrame' num2str(f)];
    print(fname,'-djpeg','-r200');
    I=imread([fname '.jpg']);
    frame=im2frame(I);
    writeVideo(M,frame);
end
end
% Sub functions
function [Dist]=mnl_MeasureEuclideanDistance(Pos1,Pos2)
Diff=Pos1-Pos2;
SqDiff=Diff.^2;
SS=sum(SqDiff);
Dist=sqrt(SS);
end
function [DistMatrix]=mnl_GroupEuD(Pos,Matrix)
%Function to measure distances to point
szM=size(Matrix);
for i=1:szM(1)
    if isnan(Matrix(i,:))==0
        [Dist]=mnl_MeasureEuclideanDistance(Pos,Matrix(i,:));
        DistMatrix(i,1)=Dist;
    else
        DistMatrix(i,1)=NaN;
    end
end
end
function [DistMatrix]=mnl_GroupEuD_v2(Pos,Matrix)
% New Method
%Calculate the distances
DistMatrix=sqrt(sum((Matrix-Pos).^2,2));
end

