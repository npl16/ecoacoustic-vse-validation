%% About
% Script to extract and plot VGGish feature embeddings from audio captured
% with Mic 1 of the 6-microphone Multichannel Autonomous Recording 
% Unit ('MAARU').

% V2.0, 03.12.2023

% This script is sctructured as follows:

% Section 1. Load data and set variables/parameters
% Section 2. Pre-process audio data and extract VGGish embedding
% Section 3. [Not featured in manuscript] Plot all 128 VGGish features over
%               time for the 4 recordings (field and 3 VSE recordings) of
%               the site specified in Section 1.
% Section 4. Plotting function 


%--------------------------------------------------------------------------
%% 1. Load VGGish CNN and all required audio, and set various parameters

net = vggish; % Load VGGissh network

numSites = 6; % Number of field sites recorded at to for which to compare field vs lab recordings
numSiteRecs = 4; % Number of times a particular site was recorded, including the original recording and lab re-recs
recToCompare = 1; % Specify which of the 4 recordings (in order: Field, Vert ReRec, 45deg ReRec, Horiz ReRec) to compare others to
numComparisons = 3; % Number of spectral differences to compute for each site (typically 3 for vertical field recording vs 3 lab re-recs in different orientations)

% *** IMPORTANT ***
sitesToStudy = 1:6; % If only studying one site in particular
siteToPlot = 3; % Specify particular site to plot in Section 3 - must be a subset of sitesToStudy

original_lengths = cell(numSites,numSiteRecs); % Cell in which to store the full length (uncropped) multi-channel audio signals (later cropped to match lengths)
fs_cell = cell(numSites,numSiteRecs); % Cell array in which to store the sampling frequencies of imported audio signals (should all be the same)
mono_sigs = cell(numSites,numSiteRecs); % Cell array in which to store cropped mono audio files from mic/channel of interest

numChans = 6; % number of channels in the spatial audio recording
ChanToUse = 1; % the specific channel being analysed in this script

nameExt = 'LP.wav';
prefixField = 'SilwoodFieldRecToUse-MC';
prefixV = 'SilwoodReRecVertToUse-MC';
prefixM = 'SilwoodReRec45ToUse-MC';
prefixH = 'SilwoodReRecHToUse-MC';
prefixes = {prefixField,prefixV,prefixM,prefixH};

orientations = {'FieldVertical','LabVertical','Lab45deg','LabHorizontal'};

% Generate filenames from prefix + orientations + extension char arrays
% above:
fileNames = cell(numSites,length(prefixes)); 

for j=1:length(prefixes)
    for i=1:numSites
        fileNames(i,j) = strcat(prefixes(j),int2str(i),nameExt);
    end
end


namesToAnalyse = fileNames(sitesToStudy,:); % Specify subset of files to import and analyse, if a subset is being used (otherwise set to full fileNames cell)
[p,q] = size(namesToAnalyse); % Get the dimensions of namesToAnalyse - corresponding to number of sites and recordings to analyse


for i=1:p
    for j=1:q
        [original_lengths{i,j}, fs_cell{i,j}] = audioread(namesToAnalyse{i,j}); %Import the audio and store in 'original_lengths' cell
    end
end

original_lengths = original_lengths(1:p,1:q); % crop 'original_lengths' to size of subset of files (specified in 'namesToAnalyse' we are working with
minLengths = min(cellfun('length',original_lengths),[],2); % Returns an array of the length of the shortest item in each row of 'original_lengths' (which contains the imported audio)

fs = fs_cell{1,1}; % use sampling frequency of first imported signal as general sampling frequency


% start and end time of signal in samples
startTime=0*fs + 1; 
idealEndTime=600*fs; 


% Crop the soundfiles for each site to match the length of the shortest recording from each site.
for i=1:p
    if minLengths(i) >= idealEndTime
        endTime = idealEndTime; % if shortest signal length is longer than idealEndTime, use idealEndTime as endTime to which to crop signal
    elseif minLengths(i) < idealEndTime
        endTime = fs*floor(minLengths(i)/fs); % if shortest signal length is shorter than idealEndTime, set endTime (in samples) to nearest second below length of shortest signal, and later crop all signals to this length.
    end
    for j=1:size(namesToAnalyse,2)
        mono_sigs{i,j} = original_lengths{i,j}(startTime:endTime,ChanToUse);
    end
end


%--------------------------------------------------------------------------
%% 2. Embed data into VGGish
% First Pre-Process Data by converting to mel-spectrograms, then pass to 
% VGGish CNN and save output in cell and as CSV

overlapPercentage = 50;
frameHopSize = 1-(overlapPercentage/100); % 'Hop' (i.e. (1-(OverlapPercent/100))*SegmentLength) between mel spectrogram 'frames' for which the 128 feature embeddings are calculated.

initWindowHop_ms = 10; % 'Hop' (i.e. (1-(OverlapPercent/100))*SegmentLength) 
% between the first set of windows/segments that the audio data is split 
% into during VGGish preprocessing. The segments are 25ms long with a 10ms 
% hop (and thus 15ms overlap).

windPerFrame = 96; % Number of these initial windows/segments that are 
% included in each mel spectrogram 'frame' for which the 128 feature 
% embeddings are calculated.

frameHopSize_ms = initWindowHop_ms*floor(windPerFrame*frameHopSize);

frameDur = 0.975; % duration in s of each mel spectrogram 'frame' 
% calculated by vggishPreprocess. The frames format the audio data in the 
% appropriate way for the vggish network to then read them.

embeddings = cell(p,q); % Cell for containing VGGish embeddings
adj_embeddings = cell(p,q); % Cell to contain VGGish embeddings for plotting
mel_specs = cell(p,q); % Cell to contain mel spectrogram frames to pass to VGGish
numFrames = cell(p,q); 
numFeatures = cell(p,q);

for i=1:p
    for j=1:q
        mel_specs{i,j} = vggishPreprocess(mono_sigs{i,j},fs,'OverlapPercentage',overlapPercentage); % Pre-process audio data into mel-spectrograms
        embeddings{i,j} = predict(net,mel_specs{i,j}); % Compute VGGish embeddings
        [numFrames{i,j},numFeatures{i,j}] = size(embeddings{i,j});
        adj_embeddings{i,j} = [embeddings{i,j};zeros(1,numFeatures{i,j})]; % For plotting, add an extra column of zeros so the final VGGish features' values appear in the top-down surface plot.
        csv_title = strcat('VGGishEmbeddings-Site',num2str(sitesToStudy),orientations{1,j},'LP.csv'); %
        % writematrix(embeddings{i,j},csv_title) % Write embeddings for current recording to CSV file
    end
end
        
numFrames = numFrames{1,1};
numFeatures = numFeatures{1,1};


%--------------------------------------------------------------------------
%% 3. Plot features for site specified in 'siteToPlot'

subplotSize = [2,2]; % Size of Subplot
plotPositions = 1:subplotSize(1)*subplotSize(2); % Array of numbers for each tile in subplot
numSubs = length(plotPositions); % Number of tiles in subplot

% Row and cell indices of adj_embeddings cell to be plotted in each tile 
rowsToPlot = [1,1,1,1]*siteToPlot;
colsToPlot = [1,2,3,4];

titles = cell(1,numSubs); % Cell to store titles of each individual plot/tile within subplot
titles{1,1} = 'Field Recording - Vertical';
titles{1,2} = 'Lab Recording - Vertical';
titles{1,3} = 'Lab Recording - 45deg';
titles{1,4} = 'Lab Recording - Horizontal';

spGlobalTitle = strcat('VGGish Feature Embeddings for (Re-)Recordings of Site',{' '},num2str(siteToPlot));

% Generate correct timestamps for labelling x-axis
% (NB. x-axis is initially y-axis, but graph is rotated by plotting
% function):
ylabIndices = [1:floor(numFrames/6):numFrames-floor(numFrames/6),numFrames];
YTicks = ylabIndices; % Locations for x-axis labels

timeToLastFrame = frameHopSize_ms*(numFrames-1);
ylabs = [0:frameHopSize_ms:timeToLastFrame,timeToLastFrame+frameDur*1000];
YTickLabs = ceil(ylabs(ylabIndices)/1000); % Values for x-axis labels

figure
for i=1:numSubs
    plotNum = [subplotSize,plotPositions(i)];
    plot_VGGish(adj_embeddings{rowsToPlot(i),colsToPlot(i)},plotNum,titles{1,i},YTicks,YTickLabs)
end
sgtitle(spGlobalTitle)


%--------------------------------------------------------------------------
%% Plotting Function

function [Y] = plot_VGGish(embeddings,plotNum,Title,YTicks,YTickLabs)
    if length(plotNum) == 1
        thisAxis = figure(plotNum);
    elseif length(plotNum) == 3
        thisAxis = subplot(plotNum(1),plotNum(2),plotNum(3));
    else
        disp("Error: plotNum must be either an integer specifying the figure numer or an 3 interger array specifying the number of rows and number of columns in the subplot, followed by the current figure's position in the suplot. See subplot documentation for details of the numbering system used.")
    end
    surf(embeddings, 'EdgeColor', 'none');
    view([90 -90])
    axis tight; 
    xlabel("Feature")
    ylabel("Time (s)")
    yticks(YTicks)
    yticklabels(YTickLabs)
    title(Title);
end