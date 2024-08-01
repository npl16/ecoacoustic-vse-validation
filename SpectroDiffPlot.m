%% About
% Script to plot spectrograms (and their differences) of recordings of a 
% pair of sites (specified by user) captured with one of the microphones on
% the 6-microphone Multichannel Autonomous Recording Unit ('MAARU').

%V6.0, 02.05.2024


% This script is structured as follows:

% Section 1. Import required sound files
% Section 2. Read and crop audio files
% Section 3. Calculate the RMS of recordings from each microphone.
% Section 4. Assign variables for pspectrum spectrogram calculation
% Section 5. Calculate Spectrograms and their differences
% Section 6. Fig. 3 - Spectral Differences only, for pair of Sites 
%            specified (3x2 subplot)
% Section 7. Fig 5 - Spectrogram and spectral differences for 
%            particular regions of pair of sites to study
% Section 8. Plotting function


%--------------------------------------------------------------------------
%% 1. Set variables and get names of audio files to import 

clear; close all

numSites = 6; % Number of field sites recorded at to for which to compare field vs lab recordings
numSiteRecs = 4; % Number of times a particular site was recorded, including the original recording and lab re-recs
recToCompare = 1; % Specify which of the 4 recordings (in order: Field, Vert ReRec, 45deg ReRec, Horiz ReRec) to compare others to
numComparisons = 3; % Number of spectral differences to compute for each site (typically 3, for vertical field recording vs 3 lab re-recs in different orientations)

% *** IMPORTANT ***
sitesToStudy = [2,5]; % Specify pair of Sites to plot
% *****************


original_lengths = cell(numSites,numSiteRecs); % Cell in which to store the full length (uncropped) audio signals (later cropped to match lengths) 
fs_cell = cell(numSites,numSiteRecs); % Cell array in which to store the sampling frequencies of imported audio signals (should all be the same, raise error if not)
sigs = cell(numSites,numSiteRecs); % Cell array in which to store cropped audio files
rmsVals = cell(numSites,numSiteRecs); % Cell array in which to store the RMS value of each microhone on the 6mic array for each recording 
rmsVar = zeros(numSites,numSiteRecs); % Matrix in which to store the variance of the RMS value of each microphone for each recording

numChans = 6; % number of channels in the spatial audio recording
ChanToUse = 1; % the specific channel to analyse / plot


% Generate names of files to import
nameExt = 'LP.wav';
prefixField = 'SilwoodFieldRecToUse-MC';
prefixV = 'SilwoodReRecVertToUse-MC';
prefixM = 'SilwoodReRec45ToUse-MC';
prefixH = 'SilwoodReRecHToUse-MC';
prefixes = {prefixField,prefixV,prefixM,prefixH};

fileNames = cell(numSites,length(prefixes));

for j=1:length(prefixes)
    for i=1:numSites
        fileNames(i,j) = strcat(prefixes(j),int2str(i),nameExt);
    end
end

namesToImport = fileNames; % Specify subset of audio files to import
namesToPlot = fileNames(sitesToStudy, :); % Specify subset of files to analyse and plot based on pair of sites set in sitesToStudy


% Set properties of colours mapped to spectrograms:
colour_range = jet(256); % full jet colour range for selecting subsets of this to use in colourmaps for spectrograms below
colours_6mic = colour_range(1:128,:); % specify first half of 'jet' colormap spectrograms of raw 6mic audio, to only get blue and green hues (as the dB values are only negative)
colours_diffs = colour_range(129:256,:); % second half for spectral differences (NB. all lab re-recordings are louder, so only the 'positive' part of the colour range is needed)

Colours = {colours_6mic, colours_diffs}; % Store colour ranges for {Spectros,Diffs} in cell
Clims = {[-200,0],[0,50]}; %Colour limits for {[Spectrograms],[SpecralDiffs]}
CAxis_Steps = [40,10]; % [Spectros,Diffs]
CBarTitles = {'Power (dB)','Difference in Power (dB)'}; % Titles for colourbars


%--------------------------------------------------------------------------
%% 2. Read and crop audio files

[p,q] = size(namesToPlot); % Get the dimensions of namesToAnalyse - corresponding to number of sites and recordings to analyse
[a,b] = size(namesToImport); % Get the dimensions of namesToImport

for i=1:a % For number of sites to analyse (NB. should always be 2)
    for j=1:b % For number of recordings (field and VSE) per site
        [original_lengths{i,j}, fs_cell{i,j}] = audioread(namesToImport{i,j}); % Import the audio and store in 'original_lengths' cell
    end
end
original_lengths = original_lengths(1:a,1:b); % crop 'original_lengths' to size of subset of files (specified in 'namesToAnalyse' we are working with
minLengths = min(cellfun('length',original_lengths),[],2); % Returns an array of the length of the shortest item in each row of 'original_lengths' (which contains the imported audio)

fs = fs_cell{1,1}; % use sampling frequency of first imported signals as general sampling frequency (TODO: need to verify that all fs in fs_cell are the same)

% start and end time of signal in samples
startTime=0*fs + 1; 
idealEndTime=600*fs; 


% Crop the soundfiles for each site to match the length of the shortest recording from each site.
for i=1:a
    if minLengths(i) >= idealEndTime
        endTime = idealEndTime; % if shortest signal length is longer than idealEndTime, use idealEndTime as endTime to which to crop signal
    elseif minLengths(i) < idealEndTime
        endTime = fs*floor(minLengths(i)/fs); % if shortest signal length is shorter than idealEndTime, set endTime (in samples) to nearest second below length of shortest signal, and later crop all signals to this length.
    end
    for j=1:size(namesToImport,2)
        sigs{i,j} = original_lengths{i,j}(startTime:endTime,:);
    end
end

sigsToPlot = sigs(sitesToStudy,:); % Specify subset of imported signals to plot


%--------------------------------------------------------------------------
%% 3. Calculate the RMS of the Signal from each Microphone for All Recordings 

for x  = 1:a % For every site...

    % Add a matrix of zeros for each site to rmsVals, with the rows 
    % corresponding to each microphone/channel from the array, and columns 
    % corresponding to each recording condition:
    % rmsVals{x,1} = zeros(numChans); 

    for y = 1:b % For every recording condition...

         % Calculate the RMS value from each channel/microphone and convert
         % to dB:
        rmsVals{x,y} = 20*log10(rms(sigs{x,y}));

        % Get variance of RMS values of each microphone:
        rmsVar(x,y) = var(rmsVals{x,y});
    end

end

% Get median, min and max of variance of RMS across microphones:
rmsVarStats = [median(rmsVar,'all'),min(min(rmsVar)),max(max(rmsVar))]; 


%--------------------------------------------------------------------------
%% 4. Assign variables for pspectrum spectrogram calculation

win_size = 0.5; % window size in seconds
M = fs*win_size; % window size in samples
win_overlap = 0.1; % proportion of window overlap
overlapPerC = win_overlap*100;
lk=0.7; % leakage of (kaiser) window


nfft = 1024; % number of DFT points for spectrogram (1024 is default for pspectrum)
L = M*win_overlap; % window overlap in samples


spectros = cell(p,q); % Cell to store matrices of spectrograms calculated from pspectrum
spectros_dB = cell(p,q); % Cell to store matrices of spectrograms calculated from pspectrum IN dB
freqs = cell(p,q);
times = cell(p,q);

spectroDiffs = cell(p,numComparisons); % Cell to store matrices of spectral differences


%--------------------------------------------------------------------------
%% 5. Calculate Spectrograms and their differences

for i=1:p % For number of sites to analyse (NB. should always be 2)
    for j=1:q % For number of recordings (field and VSE) per site

        % Use pspectrum wth parameters defined above to generate power
        % spectrogram
        [spectros{i,j}, freqs{i,j}, times{i,j}] = pspectrum(sigsToPlot{i,j}(:,ChanToUse),fs,"spectrogram",TimeResolution=win_size,...
                        OverlapPercent=overlapPerC,Leakage=lk);
        spectros_dB{i,j} = 10*log10(spectros{i,j}); % Convert to dB

        if j>=2 % For VSE re-recordings, calculate differences in spectrograms:
            spectroDiffs{i,j-1} = spectros_dB{i,j} - spectros_dB{i,recToCompare};
        end

    end
end


%--------------------------------------------------------------------------
%% 6. Fig. 3 - Spectral Differences only, for pair of Sites specified (3x2 subplot)

figure

for k = 1:2 % For two sites being plotted
    
    currentSite = k;
    
    if k == 1 % For first site being plotted
        
        % VSE Vertical Recording: 
        currentTitle = strcat('VSE Vertical vs Field - Site',{' '},int2str(sitesToStudy(k))); % Set title
        currentRec = 2; % Index for VSE Vertical spectral differences

        %Use custom plotting function
        plot_spectro(times{currentSite,currentRec}, freqs{currentSite,currentRec}, ...
            spectroDiffs{currentSite,currentRec-1}, [3,2,(k)], ...
            Clims{2}, Colours{2}, currentTitle, CBarTitles{2}, CAxis_Steps(2),[],'Frequency (Hz)');
        
        % VSE 45° Recording
        currentTitle = strcat('VSE 45° vs Field - Site',{' '},int2str(sitesToStudy(k)));
        currentRec = 3; % Index for VSE 45° spectral differences

        plot_spectro(times{currentSite,currentRec}, freqs{currentSite,currentRec}, ...
            spectroDiffs{currentSite,currentRec-1}, [3,2,(k+2)], ...
            Clims{2}, Colours{2}, currentTitle, CBarTitles{2}, CAxis_Steps(2),[],'Frequency (Hz)');
        
        % VSE Horizontal Recording: 
        currentTitle = strcat('VSE Horizontal vs Field - Site',{' '},int2str(sitesToStudy(k)));
        currentRec = 4; % Index for VSE Horizontal spectral differences

        plot_spectro(times{currentSite,currentRec}, freqs{currentSite,currentRec}, ...
            spectroDiffs{currentSite,currentRec-1}, [3,2,(k+4)], ...
            Clims{2}, Colours{2}, currentTitle, CBarTitles{2}, CAxis_Steps(2),'Time (s)','Frequency (Hz)');

    else

        % VSE Vertical Recording: 
        currentTitle = strcat('VSE Vertical vs Field - Site',{' '},int2str(sitesToStudy(k)));
        currentRec = 2; % Index for VSE Vertical spectral differences

        plot_spectro(times{currentSite,currentRec}, freqs{currentSite,currentRec}, ...
            spectroDiffs{currentSite,currentRec-1}, [3,2,(k)], ...
            Clims{2}, Colours{2}, currentTitle, CBarTitles{2}, CAxis_Steps(2),[],[]);
        
        % VSE 45° Recording: 
        currentTitle = strcat('VSE 45° vs Field - Site',{' '},int2str(sitesToStudy(k)));
        currentRec = 3; % Index for VSE Vertical spectral differences

        plot_spectro(times{currentSite,currentRec}, freqs{currentSite,currentRec}, ...
            spectroDiffs{currentSite,currentRec-1}, [3,2,(k+2)], ...
            Clims{2}, Colours{2}, currentTitle, CBarTitles{2}, CAxis_Steps(2),[],[]);
        
        % VSE Horizontal Recording: 
        currentTitle = strcat('VSE Horizontal vs Field - Site',{' '},int2str(sitesToStudy(k)));
        currentRec = 4; % Index for VSE Vertical spectral differences

        plot_spectro(times{currentSite,currentRec}, freqs{currentSite,currentRec}, ...
            spectroDiffs{currentSite,currentRec-1}, [3,2,(k+4)], ...
            Clims{2}, Colours{2}, currentTitle, CBarTitles{2}, CAxis_Steps(2),'Time (s)',[]);

    end

end

% Set overall title:
sgtitle(strcat('Spectral Differences - Sites',{' '},int2str(sitesToStudy(1)),{' '},'and',{' '},int2str(sitesToStudy(2))),'FontSize',16)

% Set figure's dimensions:
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 30, 35], 'PaperUnits', 'centimeters', 'PaperSize', [30, 35])


%--------------------------------------------------------------------------
%% 7. Fig 5 - Spectrogram and spectral differences for particular regions of pair of sites to study

regionsToStudy = {[360,420],[20,80]}; % Specify regions in seconds of two sites being analysed to plot here

for k = 1:2 % For two sites being plotted
    
    figure % Create a separate figure for each site
    currentFig = gcf;
    
    currentFig.Units = 'normalized';
    
    % YTicks = 0:100:500; 
    
    currentSite = k;
    
    newStartTimeIndex = find(ceil(times{k,1}) == regionsToStudy{k}(1), 1, 'last'); % Get the index for the start time of the region of the sound file to plot
    newEndTimeIndex = find(floor(times{k,1}) == regionsToStudy{k}(2), 1); % And for the end time
    timeArray = times{currentSite,1}(newStartTimeIndex:newEndTimeIndex,1);

    currentTitle = strcat('Field Recording Spectrogram');
    currentRec = 1; % For field recording
    plot_spectro(timeArray, freqs{currentSite,currentRec}, ... % Use custom plotting function again
        spectros_dB{currentSite,currentRec}(:,newStartTimeIndex:newEndTimeIndex), [4,1,4], ...
        Clims{1}, Colours{1}, currentTitle, CBarTitles{1}, CAxis_Steps(1),'Time (s)','Frequency (Hz)');

    currentTitle = strcat('VSE Vertical vs Field Recording');
    currentRec = 2; % For VSE Vertical recording
    plot_spectro(timeArray, freqs{currentSite,currentRec}, ...
        spectroDiffs{currentSite,currentRec-1}(:,newStartTimeIndex:newEndTimeIndex), [4,1,3], ...
        Clims{2}, Colours{2}, currentTitle, CBarTitles{2}, CAxis_Steps(2),'Time (s)','Frequency (Hz)');
    
    currentTitle = strcat('VSE 45° vs Field Recording');
    currentRec = 3; % For VSE 45° recording
    plot_spectro(timeArray, freqs{currentSite,currentRec}, ...
        spectroDiffs{currentSite,currentRec-1}(:,newStartTimeIndex:newEndTimeIndex), [4,1,2], ...
        Clims{2}, Colours{2}, currentTitle, CBarTitles{2}, CAxis_Steps(2),'Time (s)','Frequency (Hz)');
    
    currentTitle = strcat('VSE Horizontal vs Field Recording');
    currentRec = 4; % For VSE Horizontal recording
    plot_spectro(timeArray, freqs{currentSite,currentRec}, ...
        spectroDiffs{currentSite,currentRec-1}(:,newStartTimeIndex:newEndTimeIndex), [4,1,1], ...
        Clims{2}, Colours{2}, currentTitle, CBarTitles{2}, CAxis_Steps(2),'Time (s)','Frequency (Hz)');
    
    % Set overall title:
    sgtitle(strcat('Field Spectrogram and VSE Spectral Differences - Site',{' '},int2str(sitesToStudy(k))),'FontSize',17)
    
    % Set figures' dimensions:
    set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 18, 33], 'PaperUnits', 'centimeters', 'PaperSize', [18, 33])
end


%--------------------------------------------------------------------------
%% 8. Plotting Function

% Plot 6mic (Re-)Recordings' Spectrogram/Spectral Differences

function [Y] = plot_spectro(x_data,y_data,c_data,plotNum,Clims,Colours,Title,CBarTitle,CAxisSteps,XLab,YLab)

    
    if length(plotNum) == 1
        thisAxis = figure(plotNum); % Plot single pane/tile figure of plotNum is just a single number
    elseif length(plotNum) == 3 % Else create a subplot as defined by numbers in plotNum
        thisAxis = subplot(plotNum(1),plotNum(2),plotNum(3));
    else
        disp("Error: plotNum must be either an integer specifying" + ...
            " the figure numer or an 3 interger array specifying the" + ...
            " number of rows and number of columns in the subplot, " + ...
            "followed by the current figure's position in the " + ...
            "suplot. See subplot documentation for details of the " + ...
            "numbering system used.")
    end
    
    % Create 3D surface plot, top down orthogonal view
    surf(x_data, y_data, c_data, 'EdgeColor', 'none');
    axis xy;
    axis tight; 
    view(0,90);
    
    % Set colour map properties
    colourLims=Clims;
    clim(colourLims)
    colormap(thisAxis, Colours) %specify first half of 'jet' colormap for colorbar, to only get blue and green hues (as the dB values are only negative)
    
    % Set colourbar properties:
    if CBarTitle
        c = colorbar('Ticks',Clims(1):CAxisSteps:Clims(2));
        c.Label.String = CBarTitle;
        c.FontSize = 15;
    end
    
    % Set x and y axes' labels

    currentAxes = gca;
    currentAxes.XAxis.FontSize = 13;
    currentAxes.YAxis.FontSize = 13;

    if XLab
        xlabel(XLab, FontSize=15);
    end

    if YLab
        ylabel(YLab, FontSize=15);
    end
    
    % Set current sub-plot's title
    title(Title,FontSize=15);
end

