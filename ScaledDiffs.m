%% About
% Script to perform an adapted Bland Altman analysis looking at scaled
% like-for-like differences in VGGish features and 7 Acoustic Indices
% between 6mic field and VSE recordings. Also looks for trends between
% indices' differences and 6mic's recording orientation (pitch angle)

% V4.0, 02.05.2024

% This script is sctructured as follows:

% Section 1. Initialise variables; Import VGGish Embeddings and compute
%            their scaled differences.
% Section 2. Get average values and IQR of VGGish Scaled Differences.
% Section 3. Import Acoustic Indices and calculate their scaled differences
% Section 4. Get average values and IQR of Acoustic Indices' Scaled Differences
% Section 5. Supp. Fig. 1 - BA Plots of differences from all sites.
% Section 6. Figure 4 - BA Plots of overall differences across sites
% Section 7. Table 2 - Spearman's Rank Correlation Coeffs between re-recs' scaled diffs and orientation angle
% Section 8. Function for creating tiled figures of BA plots showing the  
%            medians and interquartile ranges of differences.


% Abbreviations:
% AI = Acoustic Index/Indices;
% BA = Bland Altman;
% Diffs = differences;
% Recs = recordings. 


%--------------------------------------------------------------------------
%% 1. Initialise variables; Import VGGish Embeddings and compute their scaled differences

numSites = 6; % Number of field sites recorded at for which to compare field vs lab recordings
numSiteRecs = 4; %N umber of times a particular site was recorded, including the original recording and lab re-recs

orientations = {'FieldVertical','LabVertical','Lab45deg','LabHorizontal'};
referenceRec = 1; % Specify which of the 4 recordings (in order: Field, Vert ReRec, 45deg ReRec, Horiz ReRec) to compare others to

embeddings = cell(numSites,numSiteRecs); % Cell in which to store imported VGGish embedding matrices. Rows are sites, Cols are orientations.
differencesVGG = cell(numSites,numSiteRecs-1);
scaled_differencesVGG = cell(numSites,numSiteRecs-1);
numFeatures = 128; % number of VGGish features
rangesVGG = zeros(numSites,numFeatures);

for i=1:numSites
    for j=1:numSiteRecs
        csv_title = strcat('VGGishEmbeddings-Site',num2str(i),orientations{1,j},'LP.csv'); % Get titles of CSV files with VGGish feature embeddings
        embeddings{i,j} = readmatrix(csv_title); % Store as matrices in 'embeddings' cell 
        if j==1
            rangesVGG(i,:) = range(embeddings{i,j}); % Get the range of the each of the 128 features for the field recordings
        elseif j>1
            differencesVGG{i,j-1} = embeddings{i,j} - embeddings{i,referenceRec}; % Get the differences in the 128 features for each frame between the lab re-recordings and the field recording.
            scaled_differencesVGG{i,j-1} = 100*differencesVGG{i,j-1}./rangesVGG(i,:); % Calculate the difference as a percentatge of the range of the reference data.

            % Following Becky Heath's analyses [1] looking at like-for-like
            % differences between recordings compressed to different 
            % levels, zero values where the range is 0 (which would lead to 
            % a NaN or Inf result):
            scaled_differencesVGG{i,j-1}(isnan(scaled_differencesVGG{i,j-1})) = 0; 
            scaled_differencesVGG{i,j-1}(isinf(scaled_differencesVGG{i,j-1})) = 0;
        end
    end
end


%--------------------------------------------------------------------------
%% 2. Get average values and IQR of VGGish Scaled Differences

% Initialise empty matrices to store VGGish differences' stats
avScaledDiffVGG = zeros(numSites,numSiteRecs);
medScaledDiff1VGG = zeros(numSites,numSiteRecs);
medScaledDiff2VGG = zeros(numSites,numSiteRecs);
iqrScaledDiff1VGG = zeros(numSites,numSiteRecs);
uqScaledDiffVGG = zeros(numSites,numSiteRecs);
lqScaledDiffVGG = zeros(numSites,numSiteRecs);

for p=1:numSites
    for q = 1:numSiteRecs-1
        avScaledDiffVGG(p,q+1) = mean(mean(scaled_differencesVGG{p,q}')); % Overall average scaled difference
        medScaledDiff1VGG(p,q+1) = median(mean(scaled_differencesVGG{p,q}')); % Median of mean of 128 features' differences
        medScaledDiff2VGG(p,q+1) = median(median(scaled_differencesVGG{p,q}')); % Median of median of 128 features' differences
        iqrScaledDiff1VGG(p,q+1) = iqr(mean(scaled_differencesVGG{p,q}')); % IQR of mean of mean of 128 features' differences
        uqScaledDiffVGG(p,q+1) = quantile(mean(scaled_differencesVGG{p,q}'),0.75); % Upper quartile of mean of 128 features' differences
        lqScaledDiffVGG(p,q+1) = quantile(mean(scaled_differencesVGG{p,q}'),0.25); % Lower quartile of mean of 128 features' differences
    end
end


%--------------------------------------------------------------------------
%% 3. Import Acoustic Indices and calculate their scaled differences
orientationsShort = {'FieldVert','LabVert','Lab45deg','LabH'};
AIwindow = 30; % Window duration in seconds on which acoustic indices are calculated - used in names of files containing extracted acoustic indices.

indices = cell(numSites,numSiteRecs); % Cell in which to store imported Acoustic Indices matrices. Rows are sites, Cols are orientations.
differencesAI = cell(numSites,numSiteRecs-1);
scaled_differencesAI = cell(numSites,numSiteRecs-1);
numIndices = 7; % number of Indices used
rangesAI = zeros(numSites,numIndices);

for i=1:numSites
    for j=1:numSiteRecs
        csv_title = strcat('Site',num2str(i),orientationsShort{1,j},num2str(AIwindow),'sAI-LP.csv'); %Get titles of CSV files with Acoustic Indices
        indices{i,j} = readmatrix(csv_title); % Store in as matrices in 'indices' cell 
        format long % Preserve all decimal places in CSVs.

        % Normalise the values for H index (median of the amplitude
        % envelope, in col 7) as these are based on the amplitdue of the
        % audio which varies considerably between (re-)recordings and thus
        % leads to abnormally large differences if not normalised:
        indices{i,j}(:,7) = indices{i,j}(:,7)/max(indices{i,j}(:,7));

        % And normalise Bioacoustic Index (col 4)
        indices{i,j}(:,4) = indices{i,j}(:,4)/max(indices{i,j}(:,4));

        if j==1
            rangesAI(i,:) = range(indices{i,j}); % Get the range of the each of the 7 indices for the field recordings
        elseif j>1
            differencesAI{i,j-1} = indices{i,j} - indices{i,referenceRec}; % Get the differences in the 7 indices for each frame between the lab re-recordings and the field recording.
            scaled_differencesAI{i,j-1} = 100*differencesAI{i,j-1}./rangesAI(i,:);

            % Again zero values for instances where the range is 0:
            scaled_differencesAI{i,j-1}(isnan(scaled_differencesAI{i,j-1})) = 0;
            scaled_differencesAI{i,j-1}(isinf(scaled_differencesAI{i,j-1})) = 0;
        end
    end
end


%--------------------------------------------------------------------------
%% 4. Get average values and IQR of Acoustic Indices' Scaled Differences

% Initialise empty matrices to store Acoustic Indices differences' stats
avScaledDiffAI = cell(numSites,numSiteRecs-1);
medScaledDiffAI = cell(numSites,numSiteRecs-1);
iqrScaledDiff1AI = cell(numSites,numSiteRecs-1);
uqScaledDiffAI = cell(numSites,numSiteRecs-1);
lqScaledDiffAI = cell(numSites,numSiteRecs-1);

% Cells to store the final difference statistics per index, to plot
aci = cell(numSites,1);
adi = cell(numSites,1);
aeev = cell(numSites,1);
bio = cell(numSites,1);
ndsi = cell(numSites,1);
H = cell(numSites,1);
M = cell(numSites,1);

% Create matrices to store final differences statistics per index per site 
% (to then be stored for each site in the cells above). The rows contain: 
% mean,median, IQR, upper quartile, lower quartile. Columns are for each 
% rec, compared to the original field rec. Manually add values of zero to
% the first column for 'differences' between field rec and itself. Then  
% cols r2 and up are differences for VSE Vert, VSE 45deg, and VSE Horiz.
numStats = 5; % Set number of rows in below matrices.
numPlots = 4; % Set numer of cols in below matrices.
aciCurrent = zeros(numStats,numPlots);
adiCurrent = zeros(numStats,numPlots);
aeevCurrent = zeros(numStats,numPlots);
bioCurrent = zeros(numStats,numPlots);
ndsiCurrent = zeros(numStats,numPlots);
HCurrent = zeros(numStats,numPlots);
MCurrent = zeros(numStats,numPlots);

for p=1:numSites
    for q = 1:numSiteRecs-1
        avScaledDiffAI{p,q} = mean(scaled_differencesAI{p,q});
        medScaledDiffAI{p,q} = median(scaled_differencesAI{p,q});
        iqrScaledDiff1AI{p,q} = iqr(scaled_differencesAI{p,q});
        uqScaledDiffAI{p,q} = quantile(scaled_differencesAI{p,q},0.75);
        lqScaledDiffAI{p,q} = quantile(scaled_differencesAI{p,q},0.25);

        aciCurrent(1,q+1) = avScaledDiffAI{p,q}(1);
        adiCurrent(1,q+1) = avScaledDiffAI{p,q}(2);
        aeevCurrent(1,q+1) = avScaledDiffAI{p,q}(3);
        bioCurrent(1,q+1) = avScaledDiffAI{p,q}(4);
        ndsiCurrent(1,q+1) = avScaledDiffAI{p,q}(5);
        HCurrent(1,q+1) = avScaledDiffAI{p,q}(6);
        MCurrent(1,q+1) = avScaledDiffAI{p,q}(7);

        aciCurrent(2,q+1) = medScaledDiffAI{p,q}(1);
        adiCurrent(2,q+1) = medScaledDiffAI{p,q}(2);
        aeevCurrent(2,q+1) = medScaledDiffAI{p,q}(3);
        bioCurrent(2,q+1) = medScaledDiffAI{p,q}(4);
        ndsiCurrent(2,q+1) = medScaledDiffAI{p,q}(5);
        HCurrent(2,q+1) = medScaledDiffAI{p,q}(6);
        MCurrent(2,q+1) = medScaledDiffAI{p,q}(7);

        aciCurrent(3,q+1) = iqrScaledDiff1AI{p,q}(1);
        adiCurrent(3,q+1) = iqrScaledDiff1AI{p,q}(2);
        aeevCurrent(3,q+1) = iqrScaledDiff1AI{p,q}(3);
        bioCurrent(3,q+1) = iqrScaledDiff1AI{p,q}(4);
        ndsiCurrent(3,q+1) = iqrScaledDiff1AI{p,q}(5);
        HCurrent(3,q+1) = iqrScaledDiff1AI{p,q}(6);
        MCurrent(3,q+1) = iqrScaledDiff1AI{p,q}(7);

        aciCurrent(4,q+1) = uqScaledDiffAI{p,q}(1);
        adiCurrent(4,q+1) = uqScaledDiffAI{p,q}(2);
        aeevCurrent(4,q+1) = uqScaledDiffAI{p,q}(3);
        bioCurrent(4,q+1) = uqScaledDiffAI{p,q}(4);
        ndsiCurrent(4,q+1) = uqScaledDiffAI{p,q}(5);
        HCurrent(4,q+1) = uqScaledDiffAI{p,q}(6);
        MCurrent(4,q+1) = uqScaledDiffAI{p,q}(7);

        aciCurrent(5,q+1) = lqScaledDiffAI{p,q}(1);
        adiCurrent(5,q+1) = lqScaledDiffAI{p,q}(2);
        aeevCurrent(5,q+1) = lqScaledDiffAI{p,q}(3);
        bioCurrent(5,q+1) = lqScaledDiffAI{p,q}(4);
        ndsiCurrent(5,q+1) = lqScaledDiffAI{p,q}(5);
        HCurrent(5,q+1) = lqScaledDiffAI{p,q}(6);
        MCurrent(5,q+1) = lqScaledDiffAI{p,q}(7);
    end
    
    % Store the final stats for each site in the cells for each of the
    % Acoustic Indices:
    aci{p,1} = aciCurrent;
    adi{p,1} = adiCurrent;
    aeev{p,1} = aeevCurrent;
    bio{p,1} = bioCurrent;
    ndsi{p,1} = ndsiCurrent;
    H{p,1} = HCurrent;
    M{p,1} = MCurrent;

end


%--------------------------------------------------------------------------
%% 5. Supporting Fig. 1 - BA Plots of differences from all sites

figsSoFar = 0; % Variable to count the number of figures plotted if modifying this script to do a different selection of plots.

% Plot all sites
indicesToPlot = 8; % Total number of indices (7 Acoustic Indices + averaged VGGish features differences)
subplotSize = [numSites,8]; % Will be used to create a 6-by-8 subplot
numSubs = subplotSize(1)*subplotSize(2); % Number of subplots
xTickLabs = {'Field','Vert','45°','Horiz'};
XData = [1:4];
XLims = [0,5];


subplotNum = zeros(subplotSize);
for i=1:numSites
    for j=1:indicesToPlot
        subplotNum(i,j) = j + (i-1)*indicesToPlot;
    end
end


%Set y-axix limits for Acoustic Indices to nearest 5 of upper or lower
%quartiles, or 5 if UQ less than 5 or -5 if LQ is more than -5:
for i=1:numSites

    if floor(min(lqScaledDiffVGG(i,:)/5))*5 < -5
        vggYLims(i,1) = floor(min(lqScaledDiffVGG(i,:)/5))*5;
    else
        vggYLims(i,1) = -5;
    end
    if ceil(max(uqScaledDiffVGG(i,:)/5))*5 > 5
        vggYLims(i,2) = ceil(max(uqScaledDiffVGG(i,:)/5))*5;
    else
        vggYLims(i,2) = 5;
    end

    if floor(min(aci{i,1}(5,:)/5))*5 < -5
        aciYLims(i,1) = floor(min(aci{i,1}(5,:)/5))*5;
    else
        aciYLims(i,1) = -5;
    end
    if ceil(max(aci{i,1}(4,:)/5))*5 > 5
        aciYLims(i,2) = ceil(max(aci{i,1}(4,:)/5))*5;
    else
        aciYLims(i,2) = 5;
    end

    if floor(min(adi{i,1}(5,:)/5))*5 < -5
        adiYLims(i,1) = floor(min(adi{i,1}(5,:)/5))*5;
    else
        adiYLims(i,1) = -5;
    end
    if ceil(max(adi{i,1}(4,:)/5))*5 > 5
        adiYLims(i,2) = ceil(max(adi{i,1}(4,:)/5))*5;
    else
        adiYLims(i,2) = 5;
    end

    if floor(min(aeev{i,1}(5,:)/5))*5 < -5
        aeevYLims(i,1) = floor(min(aeev{i,1}(5,:)/5))*5;
    else
        aeevYLims(i,1) = -5;
    end
    if ceil(max(aeev{i,1}(4,:)/5))*5 > 5
        aeevYLims(i,2) = ceil(max(aeev{i,1}(4,:)/5))*5;
    else
        aeevYLims(i,2) = 5;
    end
    
    if floor(min(bio{i,1}(5,:)/5))*5 < -5
        bioYLims(i,1) = floor(min(bio{i,1}(5,:)/5))*5;
    else
        bioYLims(i,1) = -5;
    end
    if ceil(max(bio{i,1}(4,:)/5))*5 > 5
        bioYLims(i,2) = ceil(max(bio{i,1}(4,:)/5))*5;
    else
        bioYLims(i,2) = 5;
    end

    if floor(min(ndsi{i,1}(5,:)/5))*5 < -5
        ndsiYLims(i,1) = floor(min(ndsi{i,1}(5,:)/5))*5;
    else
        ndsiYLims(i,1) = -5;
    end
    if ceil(max(ndsi{i,1}(4,:)/5))*5 > 5
        ndsiYLims(i,2) = ceil(max(ndsi{i,1}(4,:)/5))*5;
    else
        ndsiYLims(i,2) = 5;
    end

    if floor(min(H{i,1}(5,:)/5))*5 < -5
        HYLims(i,1) = floor(min(H{i,1}(5,:)/5))*5;
    else
        HYLims(i,1) = -5;
    end
    if ceil(max(H{i,1}(4,:)/5))*5 > 5
        HYLims(i,2) = ceil(max(H{i,1}(4,:)/5))*5;
    else
        HYLims(i,2) = 5;
    end

    if floor(min(M{i,1}(5,:)/5))*5 < -5
        MYLims(i,1) = floor(min(M{i,1}(5,:)/5))*5;
    else
        MYLims(i,1) = -5;
    end
    if ceil(max(M{i,1}(4,:)/5))*5 > 5
        MYLims(i,2) = ceil(max(M{i,1}(4,:)/5))*5;
    else
        MYLims(i,2) = 5;
    end
end


YLimsCell = {vggYLims,aciYLims,adiYLims,aeevYLims,bioYLims,ndsiYLims,HYLims,MYLims}; % Store y-axis limits in cell
indexNames = {'VGGish','ACI','ADI','AEve','Bio','NDSI','H','M'};

% Now create figure and populate subplot with custom function
figure;
currentFig = gcf;


% Set figure window size to fill screen:
currentFig.Units = 'normalized';
currentFig.Position = [0 0 1 1];
currentFig.OuterPosition = [0 0 1 1];
currentFig.WindowState = 'maximized';


% Plot each subplot:
for i=1:numSites % For each site
    for j=1:indicesToPlot % For each index to plot
        if j == 1

            title = strcat('Site',{' '},num2str(i),{' '},'-',{' '},indexNames{1,j});
            
            % Use custom plotting function (see end of script)
            BAdiffPlotV1(XData,medScaledDiff1VGG(i,:),[subplotSize,...
                subplotNum(i,j)],uqScaledDiffVGG(i,:),...
                lqScaledDiffVGG(i,:),YLimsCell{1,j}(i,:),XLims,[],title)
            currentAx = gca; % Get axis of current subplot
            currentAx.FontSize = 12; % Set overall font size to 12
            if i == numSites % For the bottom row of subplots
                set(currentAx,'XTickLabel',xTickLabs) % Include x-axis labels (for reference in all plots)
                currentAx.XTickLabelRotation = 45; % Rotate x-axis labels
            end

        elseif j == 2
            % Repeat for remaining 7 indices:
            title = strcat('Site',{' '},num2str(i),{' '},'-',{' '},indexNames{1,j});
            
            BAdiffPlotV1(XData,aci{i,1}(2,:),[subplotSize,...
                subplotNum(i,j)],aci{i,1}(4,:),aci{i,1}(5,:),...
                YLimsCell{1,j}(i,:),XLims,[],title)
            currentAx = gca;
            % For these indices, omit the y-labels (as they are already 
            % included for the leftmost column (VGGish features) of
            % subplots, and don't need to be repeated (makes things less
            % cramped):
            set(currentAx,'YLabel',[]) 
            currentAx.FontSize = 12;
            if i == numSites
                set(currentAx,'XTickLabel',xTickLabs)
                currentAx.XTickLabelRotation = 45;
            end
                    
        elseif j == 3
           
            title = strcat('Site',{' '},num2str(i),{' '},'-',{' '},indexNames{1,j});
            
            BAdiffPlotV1(XData,adi{i,1}(2,:),[subplotSize,...
                subplotNum(i,j)],adi{i,1}(4,:),adi{i,1}(5,:),...
                YLimsCell{1,j}(i,:),XLims,[],title)
            currentAx = gca;
            set(currentAx,'YLabel',[])
            currentAx.FontSize = 12;
            if i == numSites
                set(currentAx,'XTickLabel',xTickLabs)
                currentAx.XTickLabelRotation = 45;
            end
                    
        elseif j == 4
            
            title = strcat('Site',{' '},num2str(i),{' '},'-',{' '},indexNames{1,j});
            
            BAdiffPlotV1(XData,aeev{i,1}(2,:),[subplotSize,...
                subplotNum(i,j)],aeev{i,1}(4,:),aeev{i,1}(5,:),...
                YLimsCell{1,j}(i,:),XLims,[],title)
            currentAx = gca;
            set(currentAx,'YLabel',[])
            currentAx.FontSize = 12;
            if i == numSites
                set(currentAx,'XTickLabel',xTickLabs)
                currentAx.XTickLabelRotation = 45;
            end
                    
        elseif j == 5
            
            title = strcat('Site',{' '},num2str(i),{' '},'-',{' '},indexNames{1,j});
            
            BAdiffPlotV1(XData,bio{i,1}(2,:),[subplotSize,...
                subplotNum(i,j)],bio{i,1}(4,:),bio{i,1}(5,:),...
                YLimsCell{1,j}(i,:),XLims,[],title)
            currentAx = gca;
            set(currentAx,'YLabel',[])
            currentAx.FontSize = 12;
            if i == numSites
                set(currentAx,'XTickLabel',xTickLabs)
                currentAx.XTickLabelRotation = 45;
            end
                    
        elseif j == 6
            
            title = strcat('Site',{' '},num2str(i),{' '},'-',{' '},indexNames{1,j});
            
            BAdiffPlotV1(XData,ndsi{i,1}(2,:),[subplotSize,...
                subplotNum(i,j)],ndsi{i,1}(4,:),ndsi{i,1}(5,:),...
                YLimsCell{1,j}(i,:),XLims,[],title)
            currentAx = gca;
            set(currentAx,'YLabel',[])
            currentAx.FontSize = 12;
            if i == numSites
                set(currentAx,'XTickLabel',xTickLabs)
                currentAx.XTickLabelRotation = 45;
            end
                    
        elseif j == 7
            
            title = strcat('Site',{' '},num2str(i),{' '},'-',{' '},indexNames{1,j});
            
            BAdiffPlotV1(XData,H{i,1}(2,:),[subplotSize,...
                subplotNum(i,j)],H{i,1}(4,:),H{i,1}(5,:),...
                YLimsCell{1,j}(i,:),XLims,[],title)
            currentAx = gca;
            set(currentAx,'YLabel',[])
            currentAx.FontSize = 12;
            if i == numSites
                set(currentAx,'XTickLabel',xTickLabs)
                currentAx.XTickLabelRotation = 45;
            end
                    
        else
            
            title = strcat('Site',{' '},num2str(i),{' '},'-',{' '},indexNames{1,j});
            
            BAdiffPlotV1(XData,M{i,1}(2,:),[subplotSize,...
                subplotNum(i,j)],M{i,1}(4,:),M{i,1}(5,:),...
                YLimsCell{1,j}(i,:),XLims,[],title)
            currentAx = gca;
            set(currentAx,'YLabel',[])
            currentAx.FontSize = 12;
            if i == numSites
                set(currentAx,'XTickLabel',xTickLabs)
                currentAx.XTickLabelRotation = 45;
            end
            

        end
    end
end


%--------------------------------------------------------------------------
%% 6. Figure 4 - BA Plots of overall differences across sites

indicesToPlot = 8;
subplotSize = [2,4]; % For use in plotting function, to specify 2-by-4 subplot
numSubs = subplotSize(1)*subplotSize(2);
xTickLabs = {'Field','Vert','45°','Horiz'}; % x-axis labels
XData = [1:4]; % Dummy x-axis data
XLims = [0,5];
indexNames{1} = 'VGGish Features'; % Update indexNames to include add 'Features' to 'VGGish' for clarity (as it can fit on this figure)

numAIframes = height(scaled_differencesAI{1,1}); % Number of windows Acoustic Indices are calculated on

% Create cell to containing a matrix for each Acoutic index with a column
% for the each recording condition containing that index's values across
% all sites
mergedDiffsAI = cell(7,1);
for k = 1:7
    mergedDiffsAI{k,1} = zeros(numSites*numAIframes,numSiteRecs); % Initially fill the matrix with zeros
end


for a = 1:7 % For all acoustic indices:
    for b = 2:numSiteRecs % For all re-recording conditions
        for c = 1:numSites

            % Get the index range of the matrix in mergedDiffsAI that 
            % corresponds to the current site
            currentRange = ((c-1)*numAIframes + 1):c*numAIframes; 
            
            % Add index's values to corresponding part of matrix in 
            % mergedDiffsAI cell
            mergedDiffsAI{a,1}(currentRange,b) = scaled_differencesAI{c,b-1}(:,a); 
        end
    end
end

numVGGframes = height(scaled_differencesVGG{1,1}); % Number of windows VGGish features are calculated on

% Create matrix to store VGGish features' values from across all sites 
% (rows) for each recording condition (columns):
mergedDiffsAvVGG = zeros(numSites*numVGGframes,numSiteRecs);

for i = 2:numSiteRecs % For each re-recording
    for j = 1:numSites % For each site
        currentRange = ((j-1)*numVGGframes + 1):j*numVGGframes;
        mergedDiffsAvVGG(currentRange,i) = mean(scaled_differencesVGG{j,i-1}');
    end
end


% Create matrices to contain the overall medians and upper and lower
% quartiles of all indices, across all sites:
overallMedians = zeros(indicesToPlot,numSiteRecs);
overallUQs = zeros(indicesToPlot,numSiteRecs);
overallLQs = zeros(indicesToPlot,numSiteRecs);

for x = 1:numIndices+1
    if x == 1
        % The first row corresponds to the VGGIsh features' differences:
        overallMedians(x,:) = median(mergedDiffsAvVGG);
        overallUQs(x,:) = quantile(mergedDiffsAvVGG,0.75);
        overallLQs(x,:) = quantile(mergedDiffsAvVGG,0.25);
    else
        % The next 7 rows correspond to the Acoustic Indices:
        overallMedians(x,:) = median(mergedDiffsAI{x-1,:});
        overallUQs(x,:) = quantile(mergedDiffsAI{x-1,:},0.75);
        overallLQs(x,:) = quantile(mergedDiffsAI{x-1,:},0.25);
    end
end


% -------- Plot the overall medians and upper and lower quartiles: --------
figure
currentYLims = [0,0]; % Initialise vector to contain y-axis limits

for j=1:indicesToPlot % For each index to plot

    % Set y-axix limits for each index to nearest 5 of upper or lower
    % quartiles, or 5 if UQ less than 5 or -5 if LQ is more than -5:
    if floor(min(overallLQs(j,:)/5))*5 < -5
        currentYLims(1,1) = floor(min(overallLQs(j,:)/5))*5;
    else
        currentYLims(1,1) = -5;
    end

    if ceil(max(overallUQs(j,:)/5))*5 > 5
        currentYLims(1,2) = ceil(max(overallUQs(j,:)/5))*5;
    else
        currentYLims(1,2) = 5;
    end
    
    title = indexNames{1,j};
    
    % Make subplot using custom plotting function (see end of script)
    BAdiffPlotV1(XData,overallMedians(j,:),[subplotSize,...
        j],overallUQs(j,:),overallLQs(j,:),...
        currentYLims,XLims,xTickLabs,title)
                
end

% Set figure to custom size:
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 28, 14], 'PaperUnits', 'centimeters', 'PaperSize', [28, 14])


%--------------------------------------------------------------------------
%% 7. Table 2 - Spearman's Rank Correlation Coeffs between re-recs' scaled diffs and orientation angle

spearmanRes = zeros(2,8); % Store the correlation coeffs (first row) and p-values for Spearman correlations for each index (columns)

% Create matrix of zeros to store VGGish features from all sites and all
% recording consitions in one column, and the corresponding 6mic pitch 
% angle in the second column:
VGGforSpearman = zeros((numSiteRecs-1)*(height(embeddings{1,1})+...
    height(embeddings{2,1})+height(embeddings{3,1})+...
    height(embeddings{4,1})+height(embeddings{5,1})+...
    height(embeddings{6,1})),2);

heightCount = 0; % Set counter for tracking the row indieces of VGGforSpearman matrix when appending the VGGIsh differences

for a = 1:numSiteRecs-1 % For each re-recording
    for b = 1:numSites % For each site

        % Append mean scaled VGGIsh differences and orientation angle
        % to corresponding range of indices in VGGforSpearman:
        VGGforSpearman(heightCount+1:heightCount + height(scaled_differencesVGG{b,a}),1) = mean(scaled_differencesVGG{b,a}')';
        VGGforSpearman(heightCount+1:heightCount + height(scaled_differencesVGG{b,a}),2) = 0+(a-1)*45; % Calculate and append corresponding orientation angle
        heightCount = heightCount + height(scaled_differencesVGG{b,a});
    end
end

% Find Spearman's rank correlation between VGGish differences and 6mic
% orientation (pitch angle):
[rhoVgg,pvalVgg] = corr(VGGforSpearman(:,1),VGGforSpearman(:,2),'type','Spearman');
spearmanRes(:,1) = [rhoVgg;pvalVgg]; % Store in results matrix

% Create cell containing matrices of index differences (from all sites and
% recording conditions) against orientation for all indices (including
% VGGish):
AIforSpearman = cell(8,1);
AIforSpearman{1,1} = VGGforSpearman;

for k = 2:(numIndices+1) % For each acoustic index
    heightCount = 0; % re-set height counter
    AIforSpearman{k,1} = zeros((numSiteRecs-1)*numSites*numAIframes,2); % Create matrix of zeros for index differences and orientation angles
    for a = 1:numSiteRecs-1 % For all VSE re-recordings
        for b = 1:numSites % For all sites

            % Append scaled differences for each index to respective
            % matrices in AIforSpearman cell
            AIforSpearman{k,1}(heightCount+1:heightCount + height(scaled_differencesAI{b,a}),2) = scaled_differencesAI{b,a}(:,k-1);
            AIforSpearman{k,1}(heightCount+1:heightCount + height(scaled_differencesAI{b,a}),1) = 0+(a-1)*45; % Calculate and append corresponding orientation angle
            heightCount = heightCount + height(scaled_differencesAI{b,a});
        end
    end
    % Calculate Spearman's rank correlation for each index:
    [rhoAI,pvalAI] = corr(AIforSpearman{k,1}(:,1),AIforSpearman{k,1}(:,2),'type','Spearman');
    spearmanRes(:,k) = [rhoAI;pvalAI]; % Add results to spearmanRes matrix
end


%--------------------------------------------------------------------------
%% 8. BA Plots - Plotting Function

function [Y] = BAdiffPlotV1(x_data,y_data,plotNum,uq,lq,ylims,xlims,xTickLabs,Title)
    thisAxis = subplot(plotNum(1),plotNum(2),plotNum(3)); % Create subplot within arrangement described by plotNum
    hold on

    % Create background patch corresponding to ±5% limits of agreement:
    patch([0,5,5,0],[-5,-5,5,5],([115,215,150]/255),"FaceAlpha",0.5,"EdgeColor","none")
    
    % Add x-axis and other lines for other possible 'limits of agreement':
    yline(0, 'Color', 'k', 'LineWidth', 0.7) % Plot y=0 line to show x-axis
    yline(2.5, '--', 'Color', '#797979', 'LineWidth', 0.7) % Plot y=2.5
    yline(7.5, '--', 'Color', '#797979', 'LineWidth', 0.7) % Plot y=7.5
    yline(10, '--', 'Color', '#797979', 'LineWidth', 0.7) % Plot y=10
    yline(-2.5, '--', 'Color', '#797979', 'LineWidth', 0.7) % Plot y=-2.5
    yline(-7.5, '--', 'Color', '#797979', 'LineWidth', 0.7) % Plot y=-7.5
    yline(-10, '--', 'Color', '#797979', 'LineWidth', 0.7) % Plot y=-10


    scatter(x_data,y_data,72,"filled") % Plot medians
    xlim(xlims) % Set x-axis limits
    xticks([1,2,3,4]) % Set x-axis labels
    xtickangle(0)
    xticklabels(xTickLabs)

    % Plot errorbars based on upper and lower quartiles:
    errorbar(x_data,y_data,lq-y_data,uq-y_data,"LineStyle","none","LineWidth",1.5)
    ylim(ylims) % Set y-axis limits
    ylabel("% Difference") % Set y-axis label

    title(Title) % Set title
    hold off
end



