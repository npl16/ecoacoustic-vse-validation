%% About
% Script for evaluating performance of BirdNET/HARKBird for Objectives 1 
% and 2 in corresponding manuscript on using VSEs to test PAM technologies.
% 
% Here, we calculate percentage overlap of BirdNET/HARKBird
% classifications/localisations between VSE re-recordings with 6mic device
% in various orientations (i.e., different pitch angles) and 6mic field
% recordings, assess for significant differences in BirdNET/HARKBird 
% results between recording conditions, assess for correlations between  
% 6mic orientation and BirdNET/HARKBird resuts, and obtain BirdNET's
% precision and recall based on manual labelling of avian calls from the
% ambisonic field recordings of all sites.
%
% V5.0, 05.07.2024

% This script is sctructured as follows:

% Section 1. Initialise variables and data structures to store data.
% Section 2. Import BirdNET data and obtain detections per second.
% Section 3. Calculate % overlap between recordings for BirdNET detections
% Section 4. Import & calculate % overlap between recordings for HARKBird 
%            data.
% Section 5. Combine the HARKBird localisations with BirdNET detections to
%            give overall summary tables of classified and localised avian 
%            calls for each 6mic (re-)recording.
% Section 6. Look for significant differences in BirdNET confidence between 
%            recordings - presented in MS in Table 3.
% Section 7. Look for significant differences in HARKBird localisations  
%            between recordings - presented in MS in Table 3.
% Section 8. Look for Spearman's rank correlation between 6mic orientation 
%            and HARKBird localisations - presented in MS in Table 2.
% Section 9. Look for Spearman's rank correlation between 6mic orientation 
%            and BirdNET confidence.
% Section 10. Calculate Precision and Recall of BirdNET using
%             ornithologist's manual labels - presented in MS in Figure 6.
% Section 11. Create the tiled plot including BirdNET's Precision/Recall 
%             for Figure 6.


% Abbreviations:
% BN = BirdNET;
% dets = detections (i.e., BirdNET classifications)
% conf = BirdNET confidence; 
% HB = HARKBird; 
% locs = Localisations (estimated azimuth) of avian calls in HARKBird; 
% LP = low passed; 
% lspk = loudspeaker(s); 
% Recs = recordings. 


%--------------------------------------------------------------------------
%% 1. Initialise variables and data structures to store data

numSites = 6;
numRecs = 4; % Number of recordings (incl. 3 VSE recordings) for each site

BNresCell = cell(numSites,numRecs); % Cell to store raw BirdNET results (matrices to be imported from result CSV files).
BNfiltResCell = cell(numSites,numRecs); % Cell to store BirdNET results filtered by species-specific confidence thresholds.
BNlatestDets = zeros(numSites,numRecs); % Matrix to store the end time of the latest BirdNET detection from each (re-)recording.
numBNspecs = zeros(numSites,numRecs); % Matrix to store total (unfiltered) number of species detected by BirdNET in each recording.

BNscores = zeros(numSites,numRecs); % Store numerical score of number of detections in re-recs that are present in the field rec BirdNET detections.
numDetsAboveThresh = zeros(numSites,1); % Vector of number of detections in BirdNET results from the field recordings at each Site.

BNdetsPerSec = cell(numSites,numRecs); % Cell to store the detections per second - each item in the cell will be a further single dimensional cell of names of detected species (or zero if none was detected), the index of each name corresponds to its time stamp in seconds +1s to account for detections at 0s.
BNconfPerSec = cell(numSites,numRecs); % Cell to store the confidence of detections per second. As it is now known that there is only one detection per second in all recordings that is above the confidence threshold, this cell will simply contain vector arrays of the confidence of the detection per second.
BNdetsCount = cell(numSites,numRecs); % Cell to store arrays with the count of the number of detections per second

orientationNames = {'FieldRec','ReRecVert','ReRec45','ReRecH'}; % Cell of name strings for the different (re-)recording conditions/orientations
orientationAngles = [0,0,45,90]; % 6-mic array's pitch angle in recordings (clockwise angle about horizontal axis parallel to 6-mic array's plane).

genBNconfThresh = 0.5; % *Generic* confidence threshold for BirdNET predictions, used to filter out all detections below this confidence level when a species-speific confidnece threshold is not available (e.g. because the species was predicted in a VSE recording but not in the corresponding field recording).

% Create cell to contain Tables (1 for each site) with the mean field 
% recording confidence of each species detected in the field recordings:
specAvConfs = cell(numSites,1);


%--------------------------------------------------------------------------
%% 2A. Import BirdNET data

for i = 1:numSites
    for j = 1:numRecs
        tabName = strcat('Silwood',orientationNames{1,j},'ToUse-Mic1-',num2str(i),'.BirdNET.results.csv'); %Name of CSV to import
        currentTab = readtable(tabName); %read CSV file as a Table
        currentTab.Site = i*ones(height(currentTab),1);
        currentTab.Orientation = orientationAngles(j)*ones(height(currentTab),1); % Add 6mic's pitch angle that each recording was made at to every entry
        currentTab = sortrows(currentTab); % Ensure that entries are sorted by increasing time stamp
        BNresCell{i,j} = currentTab; %Add table with BN dat to cell of all tables

        % Convert species Common Names to categorical variables
        currentTab = convertvars(currentTab,"CommonName","categorical");

        % Get the mean confidence of BirdNET predictions for each species
        % in each field recording, and store in specAvConfs: 
        if j == 1 % For the field recordings

            specAvConfs{i} = varfun(@mean, currentTab, "InputVariables",...
                "Confidence", "GroupingVariables", "CommonName");
        end


        if height(currentTab) ~= 0
            BNlatestDets(i,j) = currentTab.(2)(height(currentTab)); %get the timestamp for the final BirdNET detection of a particular (re-)recording, for every site that has detections
        end

        % Count total (unfiltered) number of species detected in each
        % recording:
        numBNspecs(i,j) = height(categories(currentTab.CommonName(:)));

    end
end


finalDet = max(max(BNlatestDets)); %Get the end time of the latest detection across all recordings - the length that the 'detections per s' arrays created below should be so as to catch all possible BirdNET detections


%--------------------------------------------------------------------------
%% 2B. Convert BN data to cells of detections per second


sizes2plus = []; % Matrix in which to store the times and number of detections for each timestamp that has multiple BirdNET detections

for i = 1:numSites % For every site...

    genBNconfThresh = mean(BNresCell{i,1}.Confidence); % Set the *generic* confidence threshold to be the average confidence value of the field recording for that site.
    
    for j = 1:numRecs % For every (re-)recording

        BNdetsPerSec{i,j} = cell(finalDet+1,1); % Create a single-column cell whose length equals the latest detection time in s (we add 1 since MATLAB indexing starts at 1, not 0) in the corresponding part of the BNdetsPerSec cell for this (re-)recording
        BNdetsPerSec{i,j}(:,1) = {0}; % Fill the cell we created with a zero for each row/entry in the cell
        BNdetsCount{i,j} = zeros(finalDet+1,1); % Fill BNdetsCount with empty vectors of zeros, one for each second
        BNconfPerSec{i,j} = zeros(finalDet+1,1); % Fill BNconfPerSec with empty vectors of zeros, one for each second

        for k = 1:height(BNresCell{i,j}) % For every BirdNET detection in this (re-)recording
            
            usePrediction = 0; % Create variable to act as flag for whether species is above relevant confidence threshold

            currentTime = BNresCell{i,j}.Start_s_(k); % Get the start time of the current detection
            
            % If detection is above mean confidence FOR THAT SPECIES in field rec:
            if iscategory(specAvConfs{i}.CommonName,BNresCell{i,j}.CommonName(k))

                if BNresCell{i,j}.Confidence(k) >= specAvConfs{i}(ismember(specAvConfs{i}.CommonName, BNresCell{i,j}.CommonName(k)), 3).mean_Confidence(1)

                    usePrediction = 1;

                end

            % Else if the detection is above the generic confidence 
            % threshold (used as species not predicted for field recording):
            elseif BNresCell{i,j}.Confidence(k) >= genBNconfThresh

                usePrediction = 1;

            end


            if usePrediction == 1 % If BN prediction is above appropriate confidence threshold

                BNfiltResCell{i,j} = [BNfiltResCell{i,j};BNresCell{i,j}(k,:)];

                BNconfPerSec{i,j}(currentTime+1) = BNresCell{i,j}.Confidence(k); % Store confidence value in BNconfPerSec at index corresponding to time stamp of this BirdNET detection


                % If the cell entry at the starttime is still numeric (i.e.
                % 0; for timestamps with multiple BirdNET detections, if 
                % this is not the first detection at this time then the 
                % cell entry will be a cell with the existing detected 
                % species' name and needs to be handled differently)
                if isnumeric(BNdetsPerSec{i,j}{currentTime+1}) 
                
                    BNdetsPerSec{i,j}{currentTime+1} = BNresCell{i,j}.CommonName(k); % Pass detected species' name AS A CELL (hence round braces on (k)) to detsPerSec cell with corresponding timestamp (as its index).

                else % If there is already text in the cell (i.e. if there has(/have) already been detection(/s) at this timestamp)

                    BNdetsPerSec{i,j}{currentTime+1} = {BNdetsPerSec{i,j}{currentTime+1}{1,:},BNresCell{i,j}.CommonName{k}}; % Add this detection AS A CHAR ARRAY (hence curly braces on {k}) to the cell with existing detections at this time in BNdetsPerSec
                    sizes2plus = [sizes2plus;i,j,currentTime+1,size(BNdetsPerSec{i,j}{currentTime+1},2)]; % Add timestamp (+1 for indexing) and number of current detections for this timestamp with multiple detections as a new row in sizes2plus matrix
                end

            % Add the detection(s) added at their start time in the
            % detsPerSec cell to the next two rows so as to cover the 3s
            % detection duration of all BirdNET detections
            BNdetsPerSec{i,j}{currentTime+2} = BNdetsPerSec{i,j}{currentTime+1}; % Take whatever was at the start time and add to the next second...
            BNdetsPerSec{i,j}{currentTime+3} = BNdetsPerSec{i,j}{currentTime+1}; % ...And to the third second
            BNdetsCount{i,j}(currentTime+1:currentTime+3,1) = size(BNdetsPerSec{i,j}{currentTime+1},2);
            BNconfPerSec{i,j}(currentTime+2) = BNconfPerSec{i,j}(currentTime+1);
            BNconfPerSec{i,j}(currentTime+3) = BNconfPerSec{i,j}(currentTime+1);

            end

        end
    end
end



%% 3. Calculate % overlap between recordings for BirdNET detections


% Count number of detections in each (re-)rec (above confThresh and overall)
numDetsAboveThresh = zeros(numSites,numRecs); % Matrix of number of detections per second above their corresponding confidence threshold in BirdNET results from each (re-)recording.
numDetsOG = zeros(numSites,numRecs); % Matrix of total number of all detections per second in BirdNET results from each (re-)recording.

for i = 1:numSites
    for j = 1:numRecs
        numDetsOG(i,j) = height(BNresCell{i,j}); 
        
        for k = 1:(finalDet+1) % For every second (+1 for indexing)
            if iscell(BNdetsPerSec{i,j}{k}) % If one of the entries is a cell (i.e. contains detections) 
                numDetsAboveThresh(i,j) = numDetsAboveThresh(i,j) + size(BNdetsPerSec{i,j}{k},2); % Add the size of the cell (thus accounting for cells with multiple detections) to the count of the number of detections for this particular (re-)recording.
            end
        end

    end
end


% Tally up how many detections in the field recording of a site show up in
% the re-recordings
for i = 1:numSites
    for j = 1:numRecs
        for p = 1:(finalDet+1) % For every second (+1 for indexing)
            if isnumeric(BNdetsPerSec{i,1}{p}) == 0 & isnumeric(BNdetsPerSec{i,j}{p}) == 0 % If the fieldrec and current (re-)rec entries at this time in BNdetsPerSec both contain a detection (i.e. if both are non-numeric)
                
                for q = 1:size(BNdetsPerSec{i,j}{p},2) % For every detection at this time stamp in the (re-)rec being compared to the fieldrec...

                    % ... If this detection is in the detection(s) at this
                    % time in the corresponding field rec, then add 1 to
                    % the count of the score for this (re-)rec:
                    BNscores(i,j) = BNscores(i,j) + sum(strcmp(BNdetsPerSec{i,1}{p},BNdetsPerSec{i,j}{p}{q})); 
                    %NB. We use strcmp to compare each detection at the 
                    % current time stamp in the (re-)rec to the cell with 
                    % the character array(s) of the detection(s) at this 
                    % timestamp in the field rec.
                end

            end
        end
    end
end


% Get various stats of % overlap between (re-)recs
BNscoresPerC = 100*(BNscores./numDetsAboveThresh(:,1)); % Scores as a percentage of the number of detections above the confThresh in the field recs
BNscoresPerC(isnan(BNscoresPerC)) = 0; % If there were 0 detections in the field rec, convert these percentages (which will be NaN from dividing by 0 in the last step), to zero.

BNscoresPerC_NZ = []; % Matrix to store score values omitting sites for which each recording has a score of zero.

for i = 1:numSites
    if sum(BNscoresPerC(i,:)) ~= 0 % If all scores for this site are 0...
        BNscoresPerC_NZ = [BNscoresPerC_NZ; BNscoresPerC(i,:)]; % Omit this in scoresPerC_NZ results 
    end
end

BNmeanPerC = mean(BNscoresPerC_NZ); % Mean percentage overlap between of each each recording condition to the field recording. Cited in MS Section 5.1.3

BNscoresPerC_wReRecLocs = BNscoresPerC(find(sum(BNscoresPerC(:,2:4)')),:); % Filter down to sites that have overlapping detections in the VSE re-recordings 
BNmeanPerC_wReRecLocs = mean(BNscoresPerC_wReRecLocs); % Calculate average % overlap for sites that have overlapping detections in the VSE re-recordings

%--------------------------------------------------------------------------
%% 4. Import & calculate % overlap between recordings for HARKBird data

HBprefixes = {'FieldRec','ReRecVert','ReRec45','ReRecH'};
HBdata = cell(numSites,numRecs);

% Import and Intial processing of HB results
HBoffset = 30; % Offset between HB's 0 azimuth and the 6-mic's 0 azimuth
for i = 1:numSites
    for j = 1:numRecs
        HBdata{i,j} = readtable(strcat(HBprefixes{1,j},num2str(i),'-HBdets.csv')); % Import
        HBdata{i,j}.StartTime = floor(HBdata{i,j}.StartTime); % Round detections' start times DOWN to nearest second
        HBdata{i,j}.StartTime(HBdata{i,j}.StartTime<0) = 0; % Make any negative start times 0 (almost never happens, but is an occasional glitch in the HB outputs).
        HBdata{i,j}.EndTime = ceil(HBdata{i,j}.EndTime); % Round detections' end times UP to nearest second

        HBdata{i,j}.StartAzimuth = HBdata{i,j}.StartAzimuth - HBoffset; % Apply HBoffset to re-align recordings
        HBdata{i,j}.EndAzimuth = HBdata{i,j}.EndAzimuth - HBoffset;
        HBdata{i,j}.StartAzimuth(HBdata{i,j}.StartAzimuth<0) = HBdata{i,j}.StartAzimuth(HBdata{i,j}.StartAzimuth<0) + 360; % Convert to a 0-360° scale
        HBdata{i,j}.EndAzimuth(HBdata{i,j}.EndAzimuth<0) = HBdata{i,j}.EndAzimuth(HBdata{i,j}.EndAzimuth<0) + 360;

    end
end


HBlocs = cell(numSites,numRecs); % Cell to contain matrices of the HARKBird localisations per second
LocsInBN = cell(numSites,numRecs); % Create version of array with HB locs per second that is to be filtered to timestamps where there are BirdNET detections

HBdetsCount = cell(numSites,numRecs); % Cell to store vectors of equal length to the final BN det time, to count the number of HARKBird detections there are for each second

HBscoresOG = zeros(numSites,numRecs); % Number of original HARKBird detections for all (re-)recs - i.e, before being filtered down to those that overlap with BirdNET detections
HBdetsInBNfCount = cell(numSites,numRecs); % Number of HARKBird detections PER SECOND for all (re-)recs AFTER being filtered down to those that overlap with BirdNET fieldrec detections
HBscores = zeros(numSites,numRecs); % Number of HARKBird detections PER SECOND for all (re-)recs AFTER being filtered down to those that overlap with BirdNET fieldrec detections


%--------------------------------------------------------------------------
% Get the number of HB detections (or 'score') per second for each site and
% store in vectors in HBdetsCount cell: 

% First initialise the cell with vectors of zeros:
for i = 1:numSites
    for j = 1:numRecs
        HBdetsCount{i,j} = zeros(max(HBdata{i,j}.EndTime),1);
    end
end

% Then tally up score per second:
for i = 1:numSites
    for j = 1:numRecs

        for p = 1:height(HBdata{i,j})

            for q = (HBdata{i,j}.StartTime(p)+1):HBdata{i,j}.EndTime(p)
                HBdetsCount{i,j}(q,1) = HBdetsCount{i,j}(q,1) + 1; % Add 1 to *current score* of this time stamp (thus accounting for the 1s time stamps that appear multiple times in HB detections).
            end

        end
    end
end


%--------------------------------------------------------------------------
% Filer the HB detections for each (re-)rec down to only those that
% match times of detections in BN

for i = 1:numSites
    for j = 1:numRecs

        HBscoresOG(i,j) = sum(HBdetsCount{i,j}); % Score of number of seconds with a detection. NB. This 'score' ignores instances with multiple detections per second, as we now know that in our BirdNET detections that are above the confThresh, there are no timestamps with multiple detections.

        HBlocs{i,j} = ones(height(HBdetsCount{i,j}),max(HBdetsCount{i,j}))*370; % Initially fill HBlocs with dummy values of 370
        HBdetsInBNfCount{i,j} = zeros(finalDet+1,1);
        
        for p = 1:height(HBdata{i,j}) % For every entry in this (re-)rec's HB results
            for q = (HBdata{i,j}.StartTime(p)+1):HBdata{i,j}.EndTime(p) % For every second in the current entry
                
                k = 1; % Initialise a counter k
                
                while k <= HBdetsCount{i,j}(q,1) % While k is less than the number of times this second appears in HB detections
                    if HBlocs{i,j}(q,k) ~= 370 % If the Locations column for this second is not 370 at the current k (i.e. if the 'slot' in the current column is already taken by a previous detection of this column)...
                        k = k + 1; % Increase k (until we hit a column for this second that = 370 (i.e. a column for this second that is free))
                    
                    else
                        break
                    end
                end
                
                % Once we have a free 'slot' for this second, add the average
                % azimuth value here:
                HBlocs{i,j}(q,k) = mean([HBdata{i,j}.StartAzimuth(p),HBdata{i,j}.EndAzimuth(p)]);        
            end
        end
        

        % Create version of array with HB locs per second that is to be
        % filtered to timestamps where there are BirdNET detections:
        LocsInBN{i,j} = 370*ones(finalDet+1,size(HBlocs{i,j},2)); % Initially fill with dummy values of 370
        
        if height(HBlocs{i,j}) > (finalDet + 1) % If there are HB locs after the final BN detections
            LocsInBN{i,j} = HBlocs{i,j}([1:finalDet+1],:); % Crop to end time of latest BN detection
        else
            LocsInBN{i,j}([1:height(HBlocs{i,j})],:) = HBlocs{i,j}; % Else fill with HBlocs
        end

        for a = 1:(finalDet+1) % For every second
            if isnumeric(BNdetsPerSec{i,j}{a}) == 1 % If there are no BirdNET detections at this timestamp
                LocsInBN{i,j}(a,:) = 370; % Remove corresponding HB localisation (by overwriting with 370 dummy value)
            end
        end
        
        % Remove columns that are all values of 370 (as repeat 
        % detections of certain seconds will have been removed after 
        % filtering to just detections that match those in BN, so some 
        % extra columns may be superfluous):
        LocsInBN{i,j} = LocsInBN{i,j}(:,find(sum(abs(LocsInBN{i,j})-370))); 

    end
end


%--------------------------------------------------------------------------
% Now get score of number of HB detections from field recs that appear in 
% re-recs:

for i = 1:numSites
    for j = 1:numRecs

        % Create a copy of the field rec LocsInBN that re-recs are 
        % being compared to, so that this copied array can be 
        % manipulated as needed during the comparison for each re-rec 
        % (without modifying the 'source' reference):
        refLocs = LocsInBN{i,1}; 

        if isempty(refLocs) == 0 & isempty(LocsInBN{i,j}) == 0 % If there are detections for this fieldrec and its re-recs

            for k = 1:(finalDet+1) % For every second
                if refLocs(k,1) ~= 370 & LocsInBN{i,j}(k,1) ~= 370 % If the fieldrec and current (re-)rec entries at this time in BNdetsPerSec both contain a detection (i.e. if both are not 370)

                    p = 1; % Set counter p
                    
                    while p <= size(LocsInBN{i,j}(k,:),2) % For each detection at this second in the (re-)rec being compared:
                        q = 1; % Set / re-set value of q counter
    
                        while q <= size(refLocs(k,:),2) % For each detection at this second in the reference field rec:

                            diff = [];

                            % If there's a detection in the (re-)rec that is
                            % within ±45deg of the field rec at this time:
                            if LocsInBN{i,j}(k,p) ~= 370 && refLocs(k,q) ~= 370
                                
                                diff = LocsInBN{i,j}(k,p) - refLocs(k,q); % Get difference between current HB det value and reference HB det value
                                
                                %Then find the smallest difference between
                                %the two values
                                if abs(diff) <= 180
                                    diff = abs(diff);
                                elseif diff < -180
                                        diff = diff + 360;
                                else
                                    diff = abs(diff-360);
                                end
                            end

                            if diff <= 45 % If the difference is less than 45°
                            
                                HBdetsInBNfCount{i,j}(k,1) = HBdetsInBNfCount{i,j}(k,1) + 1; % Add one to the count of HB detections for this second of this rec that appear in the BN field dets
                                refLocs(k,q) = 370; % Set the current 'slot' (i.e. q-th detection at this second) in the reference (field rec) HB detections to 370, so that it can't be used again
                                
                                break % Break out while loop now that we have found a match
                        
                            else
                                q = q+1; % Otherwise increment q to look at next reference detection for this second
                            end

                        end
                        p = p+1; % Then increment p to look at next potential detection at this second in the (re)-rec
                    end
    
                end
            end
            
            % Get the overall score for this (re-)rec by summing the number of
            % mathcing detections per second
            HBscores(i,j) = sum(HBdetsInBNfCount{i,j}); 

        else % If the field rec or re-rec is empty (i.e. has no HB detections)
            HBscores(i,j) = 0;
        end

    end
end

HBscoresPerC = 100*(HBscores./HBscores(:,1)); % Calculate number of overalapping HB locs as % of the field locs
HBscoresPerC(isnan(HBscoresPerC)) = 0; % Zero NaN values that appear when there are no field recording locs
HBscoresPerC_NZ = HBscoresPerC(find(sum(HBscoresPerC')),:); % Filter to only sites with HB locs (i.e. remove those that are all 0)
HBmeanPerC = mean(HBscoresPerC_NZ); % Calculate average % overlap for sites with locs - cited in MS Section 5.1.3

HBscoresPerC_wReRecLocs = HBscoresPerC(find(sum(HBscoresPerC(:,2:4)')),:); % Filter down to sites that have overlapping locs in the VSE re-recordings 
HBmeanPerC_wReRecLocs = mean(HBscoresPerC_wReRecLocs); % Calculate average % overlap for sites that have overlapping locs in the VSE re-recordings


%--------------------------------------------------------------------------
%% 5. Combine the HARKBird localisations with BirdNET detections 

% 5A. First average the HB localisations to the average detection per second
% (with outliers removed)

LocsInBN2 = LocsInBN; % Create a copy of LocsInBN
avLocs = cell(numSites,numRecs); % Store the averaged LocsInBN at each second

for i = 1:numSites
    for j = 1:numRecs

        % Create array of Locs Per Sec even for recs with no filtered HB locs
        if isempty(LocsInBN2{i,j})
            LocsInBN2{i,j} = 370*ones(finalDet+1,1);
        end

        currentAvLocs = 370*ones(finalDet+1,1); % Store array of averaged Locs for each second for current rec here

        for a = 1:(finalDet+1) % For each second

            if LocsInBN2{i,j}(a,1) ~= 370 % If there are locs at this second

                % First get the magnitude of the azimuths
                LocsToCheck = LocsInBN2{i,j}(a,find(LocsInBN2{i,j}(a,:) < 370)); % Only consider non-370 values (i.e. actual HB locs)
                
                azimuthMags = abs(abs(LocsToCheck-180)-180);
                % If the most common magnitute is less than 10°, shift the values to be on a
                % scale from -180 to 180 (keeping 0 at the same point as the
                % 0-360 scale used before)
                if mode(azimuthMags) <= 10
                   LocsToCheck(find(LocsToCheck > 180)) = LocsToCheck(find(LocsToCheck > 180)) - 360;
                end
                
                % Now take average of azimuth (without outliers) and save to
                % averageLocs
                
                currentAvLocs(a,1) = mean(rmoutliers(LocsToCheck)); % Take mean of localisations with outliers removed
            end
        end

        % Now store the averaged Locs in avLocs cell
        currentAvLocs(currentAvLocs(:,1) < 0, 1) = currentAvLocs(currentAvLocs(:,1) < 0, 1) + 360; % Convert any negative values back to a 0-360° scale
        avLocs{i,j} = currentAvLocs;
    end
end



%--------------------------------------------------------------------------
% 5B. Put all BN abd HB results into single tables
% (First into tables for each site, then into overall summary tables)
%
% Process identical for all sites except site 6 due to edge cases there
% (see code below).

% Create a dummy row for adding additional rows if needed later:
dummyRowData = [0,0,0,0,370,0,0,370,0,0,370,0,0,370];

% Create table for Site 1
s1ResTab = table([0:finalDet]',1*ones(finalDet+1,1),BNdetsPerSec{1,1},BNconfPerSec{1,1},...
    avLocs{1,1},BNdetsPerSec{1,2},BNconfPerSec{1,2},avLocs{1,2},...
    BNdetsPerSec{1,3},BNconfPerSec{1,3},avLocs{1,3},BNdetsPerSec{1,4},...
    BNconfPerSec{1,4},avLocs{1,4},'VariableNames',["Time","Site","FieldDets",...
    "FieldConfs","FieldLocs","VertDets","VertConfs","VertLocs","SlantDets",...
    "SlantConfs","SlantLocs","HorizDets","HorizConfs","HorizLocs"]);
    % Note that here "Slant" refers to the 45° VSE re-recordings
    % Also NB. BN detections are in cols 3, 6, 9, and 12

% Identify rows with no detections in all recordings for deletion
toDelete = [];
for k = 1:finalDet+1
    if isnumeric(s1ResTab.(3){k}) == 1 & isnumeric(s1ResTab.(6){k}) == 1 & isnumeric(s1ResTab.(9){k}) == 1 & isnumeric(s1ResTab.(12){k}) == 1
        toDelete = [toDelete, k];
    end
end

s1ResTab(toDelete,:) = []; % Delete rows with no detections in all recordings

% Label all instances where there's no BN classification as 'None', and
% convert classifications to char arrays
for x = 1:numRecs
    for y = 1:height(s1ResTab)
        if isnumeric(s1ResTab.(3*x){y}) % If there is no BN detection
            s1ResTab.(3*x){y} = 'None';
        else
            s1ResTab.(3*x)(y) = s1ResTab.(3*x){y}; % Else convert detection to char array
        end
    end
end
s1ResTab = convertvars(s1ResTab,["FieldDets","VertDets","SlantDets","HorizDets"],"categorical");

% Table filtered to just time stamps with BN detections in the field
% recording:
s1OnlyFieldDet = s1ResTab(s1ResTab.FieldDets ~= 'None',:); 

% Further table filtered to just time stamps with matching BirdNET
% detections across all recordings and with HARKBird localisations in all
% recordings:
s1OnlyFieldDetAndLoc = s1OnlyFieldDet(s1OnlyFieldDet.FieldLocs ~= 370, :);
s1OnlyFieldDetAndLoc = s1OnlyFieldDetAndLoc(s1OnlyFieldDetAndLoc.VertLocs ~= 370, :);
s1OnlyFieldDetAndLoc = s1OnlyFieldDetAndLoc(s1OnlyFieldDetAndLoc.SlantLocs ~= 370, :);
s1OnlyFieldDetAndLoc = s1OnlyFieldDetAndLoc(s1OnlyFieldDetAndLoc.HorizLocs ~= 370, :);


% Now create table for Site 2
currentSite = 2;
s2ResTab = table([0:finalDet]',currentSite*ones(finalDet+1,1),BNdetsPerSec{currentSite,1},BNconfPerSec{currentSite,1},...
    avLocs{currentSite,1},BNdetsPerSec{currentSite,2},BNconfPerSec{currentSite,2},avLocs{currentSite,2},...
    BNdetsPerSec{currentSite,3},BNconfPerSec{currentSite,3},avLocs{currentSite,3},BNdetsPerSec{currentSite,4},...
    BNconfPerSec{currentSite,4},avLocs{currentSite,4},'VariableNames',["Time","Site","FieldDets",...
    "FieldConfs","FieldLocs","VertDets","VertConfs","VertLocs","SlantDets",...
    "SlantConfs","SlantLocs","HorizDets","HorizConfs","HorizLocs"]);

% Delete rows with no detections in all recordings
toDelete = [];
for k = 1:finalDet+1
    if isnumeric(s2ResTab.(3){k}) == 1 & isnumeric(s2ResTab.(6){k}) == 1 & isnumeric(s2ResTab.(9){k}) == 1 & isnumeric(s2ResTab.(12){k}) == 1
        toDelete = [toDelete, k];
    end
end

s2ResTab(toDelete,:) = [];

% Label all instances where there's no BN classification as 'None', and
% convert classifications to char arrays
for x = 1:numRecs
    for y = 1:height(s2ResTab)
        if isnumeric(s2ResTab.(3*x){y})
            s2ResTab.(3*x){y} = 'None';
        else
            s2ResTab.(3*x)(y) = s2ResTab.(3*x){y};
        end
    end
end
s2ResTab = convertvars(s2ResTab,["FieldDets","VertDets","SlantDets","HorizDets"],"categorical");

s2OnlyFieldDet = s2ResTab(s2ResTab.FieldDets ~= 'None',:);

s2OnlyFieldDetAndLoc = s2OnlyFieldDet(s2OnlyFieldDet.FieldLocs ~= 370, :);
s2OnlyFieldDetAndLoc = s2OnlyFieldDetAndLoc(s2OnlyFieldDetAndLoc.VertLocs ~= 370, :);
s2OnlyFieldDetAndLoc = s2OnlyFieldDetAndLoc(s2OnlyFieldDetAndLoc.SlantLocs ~= 370, :);
s2OnlyFieldDetAndLoc = s2OnlyFieldDetAndLoc(s2OnlyFieldDetAndLoc.HorizLocs ~= 370, :);


% Now create table for Site 3
currentSite = 3;
s3ResTab = table([0:finalDet]',currentSite*ones(finalDet+1,1),BNdetsPerSec{currentSite,1},BNconfPerSec{currentSite,1},...
    avLocs{currentSite,1},BNdetsPerSec{currentSite,2},BNconfPerSec{currentSite,2},avLocs{currentSite,2},...
    BNdetsPerSec{currentSite,3},BNconfPerSec{currentSite,3},avLocs{currentSite,3},BNdetsPerSec{currentSite,4},...
    BNconfPerSec{currentSite,4},avLocs{currentSite,4},'VariableNames',["Time","Site","FieldDets",...
    "FieldConfs","FieldLocs","VertDets","VertConfs","VertLocs","SlantDets",...
    "SlantConfs","SlantLocs","HorizDets","HorizConfs","HorizLocs"]);

% Delete rows with no detections in all recordings
toDelete = [];
for k = 1:finalDet+1
    if isnumeric(s3ResTab.(3){k}) == 1 & isnumeric(s3ResTab.(6){k}) == 1 & isnumeric(s3ResTab.(9){k}) == 1 & isnumeric(s3ResTab.(12){k}) == 1
        toDelete = [toDelete, k];
    end
end

s3ResTab(toDelete,:) = [];

% Label all instances where there's no BN classification as 'None', and
% convert classifications to char arrays
for x = 1:numRecs
    for y = 1:height(s3ResTab)
        if isnumeric(s3ResTab.(3*x){y})
            s3ResTab.(3*x){y} = 'None';
        else
            s3ResTab.(3*x)(y) = s3ResTab.(3*x){y};
        end
    end
end
s3ResTab = convertvars(s3ResTab,["FieldDets","VertDets","SlantDets","HorizDets"],"categorical");

s3OnlyFieldDet = s3ResTab(s3ResTab.FieldDets ~= 'None',:);

s3OnlyFieldDetAndLoc = s3OnlyFieldDet(s3OnlyFieldDet.FieldLocs ~= 370, :);
s3OnlyFieldDetAndLoc = s3OnlyFieldDetAndLoc(s3OnlyFieldDetAndLoc.VertLocs ~= 370, :);
s3OnlyFieldDetAndLoc = s3OnlyFieldDetAndLoc(s3OnlyFieldDetAndLoc.SlantLocs ~= 370, :);
s3OnlyFieldDetAndLoc = s3OnlyFieldDetAndLoc(s3OnlyFieldDetAndLoc.HorizLocs ~= 370, :);


% Now create table for Site 4
currentSite = 4;
s4ResTab = table([0:finalDet]',currentSite*ones(finalDet+1,1),BNdetsPerSec{currentSite,1},BNconfPerSec{currentSite,1},...
    avLocs{currentSite,1},BNdetsPerSec{currentSite,2},BNconfPerSec{currentSite,2},avLocs{currentSite,2},...
    BNdetsPerSec{currentSite,3},BNconfPerSec{currentSite,3},avLocs{currentSite,3},BNdetsPerSec{currentSite,4},...
    BNconfPerSec{currentSite,4},avLocs{currentSite,4},'VariableNames',["Time","Site","FieldDets",...
    "FieldConfs","FieldLocs","VertDets","VertConfs","VertLocs","SlantDets",...
    "SlantConfs","SlantLocs","HorizDets","HorizConfs","HorizLocs"]);

% Delete rows with no detections in all recordings
toDelete = [];
for k = 1:finalDet+1
    if isnumeric(s4ResTab.(3){k}) == 1 & isnumeric(s4ResTab.(6){k}) == 1 & isnumeric(s4ResTab.(9){k}) == 1 & isnumeric(s4ResTab.(12){k}) == 1
        toDelete = [toDelete, k];
    end
end

s4ResTab(toDelete,:) = [];

% Label all instances where there's no BN classification as 'None', and
% convert classifications to char arrays
for x = 1:numRecs
    for y = 1:height(s4ResTab)
        if isnumeric(s4ResTab.(3*x){y})
            s4ResTab.(3*x){y} = 'None';
        else
            s4ResTab.(3*x)(y) = s4ResTab.(3*x){y};
        end
    end
end
s4ResTab = convertvars(s4ResTab,["FieldDets","VertDets","SlantDets","HorizDets"],"categorical");

s4OnlyFieldDet = s4ResTab(s4ResTab.FieldDets ~= 'None',:);

s4OnlyFieldDetAndLoc = s4OnlyFieldDet(s4OnlyFieldDet.FieldLocs ~= 370, :);
s4OnlyFieldDetAndLoc = s4OnlyFieldDetAndLoc(s4OnlyFieldDetAndLoc.VertLocs ~= 370, :);
s4OnlyFieldDetAndLoc = s4OnlyFieldDetAndLoc(s4OnlyFieldDetAndLoc.SlantLocs ~= 370, :);
s4OnlyFieldDetAndLoc = s4OnlyFieldDetAndLoc(s4OnlyFieldDetAndLoc.HorizLocs ~= 370, :);


% Now create table for Site 5
currentSite = 5;
s5ResTab = table([0:finalDet]',currentSite*ones(finalDet+1,1),BNdetsPerSec{currentSite,1},BNconfPerSec{currentSite,1},...
    avLocs{currentSite,1},BNdetsPerSec{currentSite,2},BNconfPerSec{currentSite,2},avLocs{currentSite,2},...
    BNdetsPerSec{currentSite,3},BNconfPerSec{currentSite,3},avLocs{currentSite,3},BNdetsPerSec{currentSite,4},...
    BNconfPerSec{currentSite,4},avLocs{currentSite,4},'VariableNames',["Time","Site","FieldDets",...
    "FieldConfs","FieldLocs","VertDets","VertConfs","VertLocs","SlantDets",...
    "SlantConfs","SlantLocs","HorizDets","HorizConfs","HorizLocs"]);

% Delete rows with no detections in all recordings
toDelete = [];
for k = 1:finalDet+1
    if isnumeric(s5ResTab.(3){k}) == 1 & isnumeric(s5ResTab.(6){k}) == 1 & isnumeric(s5ResTab.(9){k}) == 1 & isnumeric(s5ResTab.(12){k}) == 1
        toDelete = [toDelete, k];
    end
end

s5ResTab(toDelete,:) = [];

% Label all instances where there's no BN classification as 'None', and
% convert classifications to char arrays
for x = 1:numRecs
    for y = 1:height(s5ResTab)
        if isnumeric(s5ResTab.(3*x){y})
            s5ResTab.(3*x){y} = 'None';
        else
            s5ResTab.(3*x)(y) = s5ResTab.(3*x){y};
        end
    end
end
s5ResTab = convertvars(s5ResTab,["FieldDets","VertDets","SlantDets","HorizDets"],"categorical");

s5OnlyFieldDet = s5ResTab(s5ResTab.FieldDets ~= 'None',:);

s5OnlyFieldDetAndLoc = s5OnlyFieldDet(s5OnlyFieldDet.FieldLocs ~= 370, :);
s5OnlyFieldDetAndLoc = s5OnlyFieldDetAndLoc(s5OnlyFieldDetAndLoc.VertLocs ~= 370, :);
s5OnlyFieldDetAndLoc = s5OnlyFieldDetAndLoc(s5OnlyFieldDetAndLoc.SlantLocs ~= 370, :);
s5OnlyFieldDetAndLoc = s5OnlyFieldDetAndLoc(s5OnlyFieldDetAndLoc.HorizLocs ~= 370, :);


%% Now create table for Site 6
% Includes additional code to handle instances of multiple BN predictions...
% at one time-stamp (only happens for Site 6 with species-specific...
% thresholding approach):

currentSite = 6;
s6ResTab = table([0:finalDet]',currentSite*ones(finalDet+1,1),BNdetsPerSec{currentSite,1},BNconfPerSec{currentSite,1},...
    avLocs{currentSite,1},BNdetsPerSec{currentSite,2},BNconfPerSec{currentSite,2},avLocs{currentSite,2},...
    BNdetsPerSec{currentSite,3},BNconfPerSec{currentSite,3},avLocs{currentSite,3},BNdetsPerSec{currentSite,4},...
    BNconfPerSec{currentSite,4},avLocs{currentSite,4},'VariableNames',["Time","Site","FieldDets",...
    "FieldConfs","FieldLocs","VertDets","VertConfs","VertLocs","SlantDets",...
    "SlantConfs","SlantLocs","HorizDets","HorizConfs","HorizLocs"]);

% Delete rows with no detections in all recordings
toDelete = [];
for k = 1:finalDet+1
    if isnumeric(s6ResTab.(3){k}) == 1 & isnumeric(s6ResTab.(6){k}) == 1 & isnumeric(s6ResTab.(9){k}) == 1 & isnumeric(s6ResTab.(12){k}) == 1
        toDelete = [toDelete, k];
    end
end

s6ResTab(toDelete,:) = [];

rowsToAdjust = cell(numRecs,1);

% Identify rows that need to be adjusted due to multiple predictions
for x = 1:numRecs
    rowsToAdjust{x,1} = [];
    for y = 1:height(s6ResTab)
        if iscell(s6ResTab.(3*x){y}) & length(s6ResTab.(3*x){y}) > 1
            if isempty(rowsToAdjust{x,1})
                rowsToAdjust{x,1} = [rowsToAdjust{x,1},y];
            elseif y > rowsToAdjust{x,1}(end)+2
                rowsToAdjust{x,1} = [rowsToAdjust{x,1},y];
            end
        end
    end
    if isempty(rowsToAdjust{x,1}) == 0
        for z = 1:length(rowsToAdjust{x,1})
            y = rowsToAdjust{x,1}(z);

            numExtraPred = length(s6ResTab.(3*x){y})-1; % Number of additional predictions
            numExtraRows = numExtraPred*3;
            BNpredToAdd = BNresCell{currentSite,x}(BNresCell{currentSite,x}.Start_s_ == s6ResTab.Time(y),:);
            BNpredToAdd = convertvars(BNpredToAdd, "CommonName", "categorical");

            % Create dummy rows in which to store the data from the
            % additional predictions:
            dummyRows = array2table(ones(numExtraRows,14).*dummyRowData);
            dummyRows.Properties.VariableNames = ["Time","Site","FieldDets",...
                "FieldConfs","FieldLocs","VertDets","VertConfs","VertLocs","SlantDets",...
                "SlantConfs","SlantLocs","HorizDets","HorizConfs","HorizLocs"];
            dummyRows.FieldDets = repmat({0},numExtraRows,1);
            dummyRows.VertDets = repmat({0},numExtraRows,1);
            dummyRows.SlantDets = repmat({0},numExtraRows,1);
            dummyRows.HorizDets = repmat({0},numExtraRows,1);
            dummyRows.Time = repmat(s6ResTab.Time(y:y+2),numExtraPred,1);
            dummyRows.Site = repmat(currentSite,numExtraRows,1);
            dummyRows.(3*x+2) = repmat(s6ResTab.(3*x+2)(y),numExtraRows,1);
            
            % Fill the dummy rows for each additional prediction
            for a = 1:numExtraPred
                currentName = s6ResTab.(3*x){y}(a);
                currentConf = BNpredToAdd(BNpredToAdd.CommonName == s6ResTab.(3*x){y}{a}, :).Confidence;
                dummyRows.(3*x)(a*1:a*1+2) = repmat({currentName},3,1);
                dummyRows.(3*x+1)(a*1:a*1+2) = repmat(currentConf,3,1);
            end

            % Now reassemble the table with the dummy rows added below
            % first prediction(s) for timestamp(s) with multiple
            % predictions:
            s6ResTab = [s6ResTab(1:y+2,:);dummyRows;s6ResTab(y+3:end,:)];
            s6ResTab.(3*x)(y:y+2) = repmat({s6ResTab.(3*x){y}(length(s6ResTab.(3*x){y}))},3,1);
        end
    end
end

% Label all instances where there's no BN classification as 'None', and
% convert classifications to char arrays
for x = 1:numRecs
    for y = 1:height(s6ResTab)
        if isnumeric(s6ResTab.(3*x){y})
            s6ResTab.(3*x){y} = 'None';
        else
            s6ResTab.(3*x)(y) = s6ResTab.(3*x){y};
        end
    end
end
s6ResTab = convertvars(s6ResTab,["FieldDets","VertDets","SlantDets","HorizDets"],"categorical");

s6OnlyFieldDet = s6ResTab(s6ResTab.FieldDets ~= 'None',:);

s6OnlyFieldDetAndLoc = s6OnlyFieldDet(s6OnlyFieldDet.FieldLocs ~= 370, :);
s6OnlyFieldDetAndLoc = s6OnlyFieldDetAndLoc(s6OnlyFieldDetAndLoc.VertLocs ~= 370, :);
s6OnlyFieldDetAndLoc = s6OnlyFieldDetAndLoc(s6OnlyFieldDetAndLoc.SlantLocs ~= 370, :);
s6OnlyFieldDetAndLoc = s6OnlyFieldDetAndLoc(s6OnlyFieldDetAndLoc.HorizLocs ~= 370, :);

%--------------------------------------------------------------------------
% 5C. Combine the various tables together: 
allResultsTab = [s1ResTab;s2ResTab;s3ResTab;s4ResTab;s5ResTab;s6ResTab]; 
OnlyFieldDets = [s1OnlyFieldDet;s2OnlyFieldDet;s3OnlyFieldDet;s4OnlyFieldDet;s5OnlyFieldDet;s6OnlyFieldDet];
OnlyFieldDetAndLocs = [s1OnlyFieldDetAndLoc;s2OnlyFieldDetAndLoc;s3OnlyFieldDetAndLoc;s4OnlyFieldDetAndLoc;s5OnlyFieldDetAndLoc;s6OnlyFieldDetAndLoc];

OnlyVertDets = allResultsTab(allResultsTab.VertDets ~= 'None',:);
OnlySlantDets = allResultsTab(allResultsTab.SlantDets ~= 'None',:);
OnlyHorizDets = allResultsTab(allResultsTab.HorizDets ~= 'None',:);


%--------------------------------------------------------------------------
%% 6. Table 3 - Look for significant differences in BN confidence between recordings

% Create matrix with the confidence values of all BN detections filtered to
% those that appear in all recording, and run Friedman test with Dunn 
% post-hoc to compare between groups:

filteredAllRes = allResultsTab(allResultsTab.FieldConfs ~= 0, :);
filteredAllRes = filteredAllRes(filteredAllRes.VertConfs ~= 0, :);
filteredAllRes = filteredAllRes(filteredAllRes.SlantConfs ~= 0, :);
filteredAllRes = filteredAllRes(filteredAllRes.HorizConfs ~= 0, :);
filteredConfMat = [filteredAllRes.FieldConfs,filteredAllRes.VertConfs,filteredAllRes.SlantConfs,filteredAllRes.HorizConfs];


% Get median and standard deviation of confidence values for predictions
% that appear across all recording conditions:
filteredConfMatMedian = median(filteredConfMat)
filteredConfMatStd = std(filteredConfMat)

% Plot histograms of the confidence values for BN detections in all recs
figure
t = tiledlayout(2,2,'TileSpacing','loose','Padding','compact');
nexttile(1)
histfit(filteredConfMat(:,1)) % Field recording confidence values
title('Field Recs')
nexttile(2)
histfit(filteredConfMat(:,2)) % VSE Vertical recording confidence values
title('VSE Vert Recs')
nexttile(3)
histfit(filteredConfMat(:,3)) % VSE 45° ('slant') recording confidence values
title('VSE 45 Recs')
nexttile(4)
histfit(filteredConfMat(:,4)) % VSE Horizontal recording confidence values
title('VSE Horiz Recs')
title(t, 'Histograms of BN Confidence')

[Pfilt,Tfilt,Sfilt] = friedman(filteredConfMat) %Friedman test
[COMPfilt,MEANSfilt,Hfilt,GNAMESfilt] = multcompare(Sfilt,'ctype','dunn-sidak') % Dunn post-hoc - data presented in Table 3


%--------------------------------------------------------------------------
%% 7. Table 3 - Look for significant differences in HB locs between recordings

% Plot histograms of HB locs for birds localised in all recs
figure
t = tiledlayout(2,2,'TileSpacing','loose','Padding','compact');
nexttile(1)
histfit(OnlyFieldDetAndLocs.FieldLocs)
title('Field Recs')
nexttile(2)
histfit(OnlyFieldDetAndLocs.VertLocs)
title('VSE Vert Recs')
nexttile(3)
histfit(OnlyFieldDetAndLocs.SlantLocs)
title('VSE 45 Recs')
nexttile(4)
histfit(OnlyFieldDetAndLocs.HorizLocs)
title('VSE Horiz Recs')
title(t, 'Histograms of HB Localisations')

% Put HB locs for birds localised in all recs into matrix
locsMat = [OnlyFieldDetAndLocs.FieldLocs, OnlyFieldDetAndLocs.VertLocs, OnlyFieldDetAndLocs.SlantLocs, OnlyFieldDetAndLocs.HorizLocs];

% Get median and standard deviation of HARKBird localisations that appear
% across all recording conditions:
locsMatMedian = median(locsMat)
locsMatStd = std(locsMat)

% As data are not normally distributed and dependent, use Friedman test
% with Dunn post-hoc again
[Plocs,Tlocs,Slocs] = friedman(locsMat)
[COMPlocs,MEANSlocs,Hlocs,GNAMESlocs] = multcompare(Slocs,'ctype','dunn-sidak') % Dunn post-hoc - data presented in Table 3


%--------------------------------------------------------------------------
%% 8. Table 2 - Look for Spearman's rank correlation between 6mic orientation and HB loc (Table 2)

% Create 2-column matrix, first column containing 6mic orientation angle
% for each HB localisation, and second column containing HB locs for birds
% detected in all recording conditions (which we can obtain from locsMat,
% created in the previous section):
HBlocsForSpearman = zeros(height(locsMat)*3,2); 
% As we will effectively be taking the columns 2, 3, and 4 of locsMat and 
% concatenating them vertically in the second column, this matrix must be
% 3x the height of locsMat

for x=1:3 % For each re-rec condition
    startPos = (x-1)*height(locsMat) + 1; % Get the start index in HBlocsForSpearman for this condition
    endPos = (x-1)*height(locsMat) + height(locsMat); % Get the end index in HBlocsForSpearman for this condition

    % Add vector of orientation angle to first column, and HB locs from 
    % that orientation to the second column
    HBlocsForSpearman(startPos:endPos,:) = [orientationAngles(x+1)*ones(height(locsMat),1),locsMat(:,x+1)]; 
end

[rhoLocs,pvalLocs] = corr(HBlocsForSpearman(:,1), HBlocsForSpearman(:,2), 'type','Spearman');


%--------------------------------------------------------------------------
%% 9. Table 2 - Look for Spearman's rank correlation between 6mic orientation and BirdNET confidence (Table 2)

% Compile BirdNET results of all VSE re-recordings into single table:
allReRecs = [BNfiltResCell{1,2} ; BNfiltResCell{2,2} ; BNfiltResCell{3,2} ; ...
    BNfiltResCell{4,2} ; BNfiltResCell{5,2} ; BNfiltResCell{6,2} ; BNfiltResCell{1,3} ; ...
    BNfiltResCell{2,3} ; BNfiltResCell{3,3} ; BNfiltResCell{4,3} ; BNfiltResCell{5,3} ; ...
    BNfiltResCell{6,3} ; BNfiltResCell{1,4} ; BNfiltResCell{2,4} ; BNfiltResCell{3,4} ; ...
    BNfiltResCell{4,4} ; BNfiltResCell{5,4} ; BNfiltResCell{6,4}];


allReRecs = convertvars(allReRecs,["ScientificName","CommonName"],"categorical"); % Convert species' names to categorical variables
detsByBird = groupcounts(allReRecs,"CommonName"); % Get the frequency of each BN detection...
detsByBird = flip(sortrows(detsByBird,'GroupCount')); % ... and order this in terms of highest to lowest frequency

% Plot frequeny distribution of BN detections of each species
figure
bar(detsByBird.CommonName,detsByBird.GroupCount)

% Calculate Spearman's rank correlation coefficient between 6mic orientation and
% the confidence values of all BN detections
[rho,pval] = corr(allReRecs.Orientation,allReRecs.Confidence,'type','Spearman');

% Store the result in a new table:
BNSpearmanResults = table("All",rho,pval); 
BNSpearmanResults.Properties.VariableNames = ["CommonName","Rho","p-Value"];
BNSpearmanResults = convertvars(BNSpearmanResults,"CommonName","categorical");

for k = 1:4 % For each of the 4 most frequently detected bird species:
    tabForSpearman = allReRecs(allReRecs.CommonName == detsByBird.CommonName(k),:); % Extract the species' BirdNET results
    [rhoCurrent,pvalCurrent] = corr(tabForSpearman.Orientation,tabForSpearman.Confidence,'type','Spearman'); % Get Spearman's rank correlation
    BNSpearmanResults = [BNSpearmanResults;{detsByBird.CommonName(k),rhoCurrent,pvalCurrent}]; % Add this to the results table, which feeds into Table 2 in the MS
end


%--------------------------------------------------------------------------
%% 10. Precision and Recall of BirdNET (Fig. 6)


BNthreshs = [0, 0.25, 0.5, 0.75, 0.9];
numThresh = length(BNthreshs);
BNresMultThresh = cell(5,1);
BNresMultThresh(:,1) = {BNresCell;BNresCell;BNresCell;BNresCell;BNresCell};

for count = 1:5
    for m = 1:6
        for n = 1:4
            BNresMultThresh{count,1}{m,n} = BNresMultThresh{count,1}{m,n}(BNresMultThresh{count,1}{m,n}.Confidence > BNthreshs(count),:);
        end
    end
end



%% Existing precision and recall script starts from here

% Initialise further variables and data structures for storing manual
% labels and more:

manualLabels = cell(numSites,1); % Cell to store tables of Manual Labels imported from CSV files
numWindows = 200; % Number of 3 s windows per recording used for BirdNET detections and Manual Labels
detsPerWindow = cell(numSites,numThresh); % Cell to store counts the number of (unfiltered) BN detections per 3 s window of each recording for each site, for each BN conf threshold (cell columns)
MLflags = cell(numSites,numThresh); % Cell to store tables with rows for each time window and cols for TP/TN/FP/FN, for each BN conf threshold (cell columns)

MLstats = cell(numSites,numThresh); % Cell to store matrices with the precision and recall of BN of each recording for each site, for each BN conf threshold (cell columns)

numManDets = zeros(numSites,1); % Array to store total number of manual labels for each recording - Presented in Figure 6
manSpecs = cell(numSites,1); % Cell to store names of manually-labelled species from each site
numManSpecs = zeros(numSites,1); % Matrix to store number of manually-labelled species in each recording  - Presented in Figure 6


% Import Ornithologist's manual species counts:
manualLabels{1,1} = readtable('SilwoodRecsS1-ManualLabelsV1Complete.csv','ReadVariableNames',true);
manualLabels{2,1} = readtable('SilwoodRecsS2-ManualLabelsV1Complete.csv','ReadVariableNames',true);
manualLabels{3,1} = readtable('SilwoodRecsS3-ManualLabelsV1Complete.csv','ReadVariableNames',true);
manualLabels{4,1} = readtable('SilwoodRecsS4-ManualLabelsV1Complete.csv','ReadVariableNames',true);
manualLabels{5,1} = readtable('SilwoodRecsS5-ManualLabelsV1Complete.csv','ReadVariableNames',true);
manualLabels{6,1} = readtable('SilwoodRecsS6-ManualLabelsV1Complete.csv','ReadVariableNames',true);
% The CSV of manaual labels contains a vector for all 3 s timestamps within
% a 10min period (max. duration of the recorings) and vectors with entries
% of species' names for time stamps where those species were heard


for s = 1:numSites % For all sites

    % If one of the Species' columns in the Manual Labels' CSVs is empty
    % all values in that column will be read as NaN during import, so check
    % first if first value in each Species name column is a cell (rather 
    % than NaN, and replace entire column with empty char array ('') if so:
    if iscell(manualLabels{s,1}.Species1(1)) == 0
        manualLabels{s,1}.Species1 = repmat({''},200,1);
    end
    if iscell(manualLabels{s,1}.Species2(1)) == 0
        manualLabels{s,1}.Species2 = repmat({''},200,1);
    end
    if iscell(manualLabels{s,1}.Species3(1)) == 0
        manualLabels{s,1}.Species3 = repmat({''},200,1); 
    end

    % Collapse instances with multiple species to all be in the leftmost free 
    % position for species names in the table
    for i = 1:numWindows

        if isempty(manualLabels{s,1}.Species3{i}) == 0 & isempty(manualLabels{s,1}.Species2{i}) == 1
            manualLabels{s,1}.Species2{i} = manualLabels{s,1}.Species3{i};
            manualLabels{s,1}.Species3{i} = '';
        end
        if isempty(manualLabels{s,1}.Species2{i}) == 0 & isempty(manualLabels{s,1}.Species1{i}) == 1
            manualLabels{s,1}.Species1{i} = manualLabels{s,1}.Species2{i};
            manualLabels{s,1}.Species2{i} = '';
        end
    end
    
    convertvars(manualLabels{s,1},["Species1","Species2","Species3"],"categorical");


    % Extract relevant columns and convert species' names to categorical
    % variables:
    manualLabels{s,1} = manualLabels{s,1}(:,[1:5]);
    
    % Put all manual labels for the current site into a single-column table
    % to then count how many species there are
    currentManLabs = table([manualLabels{s,1}.Species1(:);manualLabels{s,1}.Species2(:);manualLabels{s,1}.Species3(:)],'VariableNames',"ManualLabels");
    currentManLabs = convertvars(currentManLabs,["ManualLabels"],"categorical");

    % Store names of manually-labelled species for this site in
    % manSpecs, and number of manually-labelled species in numManSpecs:
    manSpecs{s,1} = categories(currentManLabs.ManualLabels);
    numManSpecs(s) = height(manSpecs{s,1});  % Presented in Figure 6

    % Create additional column in manualLabels containing cell array of 
    % cells with species' names (so all species at each timestamp can be 
    % accessed in one cell of the table) 
    dummyCell = cell(numWindows,1); % First make a dummy cell...
    dummyCell(:,1) = {0}; % ... and fill it with zeros.
    
    manualLabels{s,1}.Manual_Det_Array = dummyCell; % Append this to the table
    
    for i = 1:numWindows
        if isempty(manualLabels{s,1}.Species1{i}) == 0 % If there is a manual detection
            % Add first species AS A CELL (round braces on 'Species1(i)'):
            manualLabels{s,1}.Manual_Det_Array{i} = manualLabels{s,1}.Species1(i);
            numManDets(s) = numManDets(s) + 1; % Add 1 to count of number of manual labels for this site - Presented in Figure 6
        end
        
        if isempty(manualLabels{s,1}.Species2{i}) == 0 % If there is a second manual detection
            % Append second species to cell array AS A CELL (curly braces on
            % 'Species2.{i}'):
            manualLabels{s,1}.Manual_Det_Array{i} = {manualLabels{s,1}.Manual_Det_Array{i}{1,:},manualLabels{s,1}.Species2{i}};
            numManDets(s) = numManDets(s) + 1; % Add 1 to count of number of manual labels for this site - Presented in Figure 6
        end
        
        if isempty(manualLabels{s,1}.Species3{i}) == 0 % If there is a third manual detection
            % Append third species to cell array AS A CELL (curly braces on
            % 'Species3.{i}'):
            manualLabels{s,1}.Manual_Det_Array{i} = {manualLabels{s,1}.Manual_Det_Array{i}{1,:},manualLabels{s,1}.Species3{i}};
            numManDets(s) = numManDets(s) + 1; % Add 1 to count of number of manual labels for this site - Presented in Figure 6
        end
    end

end


for threshCount = 1:numThresh % For every BN confidence threshold

    for s = 1:numSites % For all sites
    
        % Count the number of (unfiltered) BN detections per 3 s window
        detsPerWindow{s,threshCount} = zeros(numWindows,numRecs); % First create matrix of zeros
        
        for a = 1:numRecs % For each recording
            for b = 1:height(BNresMultThresh{threshCount,1}{s,a}) % For each BN detection in the current recording
                currentWindow = (BNresMultThresh{threshCount,1}{s,a}.Start_s_(b)/3) + 1; % Index of current window is detection's start time divided by 3 and +1
                detsPerWindow{s,threshCount}(currentWindow,a) = detsPerWindow{s,threshCount}(currentWindow,a) + 1; % Add detection to count of detections per window
            end
        end
        
        
        % Create cell (one column for each recording condition) to store flags of 
        % whether the BN detection (or lack thereof) in each 3 s window is a 
        % True Positive, True Negative, False Positive, or False Negative:
        MLflags{s,threshCount} = cell(1,numRecs); 
        
        % Store tables with TP, TN, FP, and FN in each entry of the cell
        for k = 1:numRecs
            MLflags{s,threshCount}{1,k} = table(zeros(numWindows,1),zeros(numWindows,1),...
                zeros(numWindows,1),zeros(numWindows,1),'VariableNames',["TP",...
                "TN","FP","FN"]);
        end
        
        % Create empty matrix to store the Precision (top row) and Recall (bottom 
        % row) of BN on each 6mic recording of Site 2:
        MLstats{s,threshCount} = zeros(2,numRecs); % Presented in Figure 6
        
        
        for i = 1:numRecs % For each 6mic recording
            for j = 1:numWindows % For each 3 s window
                
                % If there is(/are) a manual detection(/s) and BN detection(/s) at
                % the current 3 s window:
                if isempty(manualLabels{s,1}.Species1{j}) == 0 & detsPerWindow{s,threshCount}(j,i) ~= 0
                                
                    for p = 1:detsPerWindow{s,threshCount}(j,i) % For the number of BN detections at this time stamp:
                        
                        % Get the index of the current windows in BNresCell (and add
                        % p-1 to account for windows with multiple detections, as 
                        % each detection, even if for the same window, is listed as
                        % a new row in BNresCell)
                        indexOfBNdet = find(BNresMultThresh{threshCount,1}{s,i}.Start_s_ == manualLabels{s,1}.StartTime_s_(j),1) + p - 1;
                        
                        % Use strcmp to compare the array of cells containing the
                        % common names of manually-detected species to the 
                        % char array with the common name of species detected by
                        % BN in this window:
                        res = strcmp(manualLabels{s,1}.Manual_Det_Array{j},BNresMultThresh{threshCount,1}{s,i}.CommonName{indexOfBNdet});
                        % This outputs an array the same length as the current
                        % Manual_Det_Array with a 1 in the position of the matching
                        % common name (if there is one) and 0s elsewhere
            
                        if sum(res) == 0 % If there is no matching detecion:
                            MLflags{s,threshCount}{1,i}.FP(j) = MLflags{s,threshCount}{1,i}.FP(j) + 1; % False Positive
                        else 
                            MLflags{s,threshCount}{1,i}.TP(j) = MLflags{s,threshCount}{1,i}.TP(j) + 1; % True Positive
                        end
                            
                    end
            
                    % If there are fewer BN dections for this timestamp than there
                    % are manual detections, log the difference between them as 
                    % False Negatives:
                    DiffInDets = length(manualLabels{s,1}.Manual_Det_Array{j}) - detsPerWindow{s,threshCount}(j,i);
                    
                    if DiffInDets > 0
                        MLflags{s,threshCount}{1,i}.FN(j) = DiffInDets; % False Negative(s)
                    end
        
                elseif isempty(manualLabels{s,1}.Species1{j}) == 0 & detsPerWindow{s,threshCount}(j,i) == 0 % If manual det but no BN det
                    
                    MLflags{s,threshCount}{1,i}.FN(j) = length(manualLabels{s,1}.Manual_Det_Array{j}); % False Negative(s) 
            
                elseif isempty(manualLabels{s,1}.Species1{j}) == 1 & detsPerWindow{s,threshCount}(j,i) ~= 0 % If no manual det but BN det
            
                    MLflags{s,threshCount}{1,i}.FP(j) = detsPerWindow{s,threshCount}(j,i); % False Positive(s)
            
                else
            
                    MLflags{s,threshCount}{1,i}.TN(j) = 1; % True Negative
                end
            
            end
        
            MLstats{s,threshCount}(1,i) = sum(MLflags{s,threshCount}{1,i}.TP) / (sum(MLflags{s,threshCount}{1,i}.TP) + sum(MLflags{s,threshCount}{1,i}.FP)); % Calculate precision - presented in Figure 6
            MLstats{s,threshCount}(2,i) = sum(MLflags{s,threshCount}{1,i}.TP) / (sum(MLflags{s,threshCount}{1,i}.TP) + sum(MLflags{s,threshCount}{1,i}.FN)); % Calculate recall - presented in Figure 6
        
        end
    
    end

end


%--------------------------------------------------------------------------
%% 11. Figure 6 - Create Plot of BirdNET's Precision and Recall

% Colour blind friendly chart colours (from Okabe and Ito) in RGB:
chartColours = [0 158 115; 0 114 178; 86 180 233; 240 228 66; 204 121 167; 230 159 0;  0 0 0; 213 94 0]/255;

% Create cell to store matrices with number of BN species predictions...
% for each recording at various confidence thresholds:
nSpecsMultThresh = cell(numThresh,1);

% And a cell for the number of calls detected by BN at each conf. thresh.:
nCallsMultThresh = cell(numThresh,1);

for a = 1:numThresh % For each confidence threshold...

    % Fill cells with matrices of zeroes: 
    nSpecsMultThresh{a,1} = zeros(numSites,numRecs); 
    nCallsMultThresh{a,1} = zeros(numSites,numRecs);

    for b = 1:numSites % For each site... 
        for c = 1:numRecs % For each (re-)rec...

            nCallsMultThresh{a,1}(b,c) = height(BNresMultThresh{a,1}{b,c}); % Get the number of calls

            % Store version of filtered BN results table with CommonNames...
            % as categorical in temporary table (orignical BNresMultThresh...
            % must not have categorical variables for earlier processing):
            currentBNres = convertvars(BNresMultThresh{a,1}{b,c}, "CommonName", "categorical"); 

            nSpecsMultThresh{a,1}(b,c) = height(categories(currentBNres.CommonName)); % Get the number of species

        end
    end
end


% Plot figure -------------------------------------------------------------
figure
tiledlayout(4,1)
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 27, 30], 'PaperUnits', 'centimeters', 'PaperSize', [27, 30]) % Set dimensions
Xlims = [0.5 4.5];
siteXpos = linspace(0.75, 1.25, 6); % relative x data positions for first plot

ax1 = nexttile(1);  % Precision and Recall --------------------------------

%Site 1
site = 1;
LWidth = 1.25;

pre1_3 = scatter(siteXpos(site)+[0:3], MLstats{site,3}(1,:), 'o', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth);
hold on
pre1_2 = scatter(siteXpos(site)+[0:3], MLstats{site,2}(1,:), 'v', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth);
pre1_4 = scatter(siteXpos(site)+[0:3], MLstats{site,4}(1,:), '^', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth);
pre1_3.DisplayName = 'Site 1 Precision';
pre1_2.DisplayName = 'Site 1 (0.25)';
pre1_4.DisplayName = 'Site 1 (0.75)';

rec1_3 = scatter(siteXpos(site)+[0:3], MLstats{site,3}(2,:), 'o', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth, "MarkerFaceColor", chartColours(site,:), "MarkerFaceAlpha", 0.3);
rec1_2 = scatter(siteXpos(site)+[0:3], MLstats{site,2}(2,:), 'v', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth, "MarkerFaceColor", chartColours(site,:), "MarkerFaceAlpha", 0.3);
rec1_4 = scatter(siteXpos(site)+[0:3], MLstats{site,4}(2,:), '^', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth, "MarkerFaceColor", chartColours(site,:), "MarkerFaceAlpha", 0.3);
rec1_3.DisplayName = 'Site 1 Recall';

% Site 2
site = 2;
pre2_3 = scatter(siteXpos(site)+[0:3], MLstats{site,3}(1,:), 'o', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth);
pre2_2 = scatter(siteXpos(site)+[0:3], MLstats{site,2}(1,:), 'v', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth);
pre2_4 = scatter(siteXpos(site)+[0:3], MLstats{site,4}(1,:), '^', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth);
pre2_3.DisplayName = 'Site 2 Precision';
pre2_2.DisplayName = 'Site 2 (0.25)';
pre2_4.DisplayName = 'Site 2 (0.75)';
rec2_3 = scatter(siteXpos(site)+[0:3], MLstats{site,3}(2,:), 'o', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth, "MarkerFaceColor", chartColours(site,:), "MarkerFaceAlpha", 0.3);
rec2_2 = scatter(siteXpos(site)+[0:3], MLstats{site,2}(2,:), 'v', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth, "MarkerFaceColor", chartColours(site,:), "MarkerFaceAlpha", 0.3);
rec2_4 = scatter(siteXpos(site)+[0:3], MLstats{site,4}(2,:), '^', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth, "MarkerFaceColor", chartColours(site,:), "MarkerFaceAlpha", 0.3);
rec2_3.DisplayName = 'Site 2 Recall';

% Site 3
site = 3;
pre3_3 = scatter(siteXpos(site)+[0:3], MLstats{site,3}(1,:), 'o', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth);
pre3_2 = scatter(siteXpos(site)+[0:3], MLstats{site,2}(1,:), 'v', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth);
pre3_4 = scatter(siteXpos(site)+[0:3], MLstats{site,4}(1,:), '^', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth);
pre3_3.DisplayName = 'Site 3 Precision';
pre3_2.DisplayName = 'Site 3 (0.25)';
pre3_4.DisplayName = 'Site 3 (0.75)';
rec3_3 = scatter(siteXpos(site)+[0:3], MLstats{site,3}(2,:), 'o', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth, "MarkerFaceColor", chartColours(site,:), "MarkerFaceAlpha", 0.3);
rec3_2 = scatter(siteXpos(site)+[0:3], MLstats{site,2}(2,:), 'v', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth, "MarkerFaceColor", chartColours(site,:), "MarkerFaceAlpha", 0.3);
rec3_4 = scatter(siteXpos(site)+[0:3], MLstats{site,4}(2,:), '^', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth, "MarkerFaceColor", chartColours(site,:), "MarkerFaceAlpha", 0.3);
rec3_3.DisplayName = 'Site 3 Recall';

% Site 4
site = 4;
pre4_3 = scatter(siteXpos(site)+[0:3], MLstats{site,3}(1,:), 'o', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth);
pre4_2 = scatter(siteXpos(site)+[0:3], MLstats{site,2}(1,:), 'v', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth);
pre4_4 = scatter(siteXpos(site)+[0:3], MLstats{site,4}(1,:), '^', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth);
pre4_3.DisplayName = 'Site 4 Precision';
pre4_2.DisplayName = 'Site 4 (0.25)';
pre4_4.DisplayName = 'Site 4 (0.75)';
rec4_3 = scatter(siteXpos(site)+[0:3], MLstats{site,3}(2,:), 'o', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth, "MarkerFaceColor", chartColours(site,:), "MarkerFaceAlpha", 0.3);
rec4_2 = scatter(siteXpos(site)+[0:3], MLstats{site,2}(2,:), 'v', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth, "MarkerFaceColor", chartColours(site,:), "MarkerFaceAlpha", 0.3);
rec4_4 = scatter(siteXpos(site)+[0:3], MLstats{site,4}(2,:), '^', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth, "MarkerFaceColor", chartColours(site,:), "MarkerFaceAlpha", 0.3);
rec4_3.DisplayName = 'Site 4 Recall';

% Site 5
site = 5;
pre5_3 = scatter(siteXpos(site)+[0:3], MLstats{site,3}(1,:), 'o', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth);
pre5_2 = scatter(siteXpos(site)+[0:3], MLstats{site,2}(1,:), 'v', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth);
pre5_4 = scatter(siteXpos(site)+[0:3], MLstats{site,4}(1,:), '^', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth);
pre5_3.DisplayName = 'Site 5 Precision';
pre5_2.DisplayName = 'Site 5 (0.25)';
pre5_4.DisplayName = 'Site 5 (0.75)';
rec5_3 = scatter(siteXpos(site)+[0:3], MLstats{site,3}(2,:), 'o', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth, "MarkerFaceColor", chartColours(site,:), "MarkerFaceAlpha", 0.3);
rec5_2 = scatter(siteXpos(site)+[0:3], MLstats{site,2}(2,:), 'v', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth, "MarkerFaceColor", chartColours(site,:), "MarkerFaceAlpha", 0.3);
rec5_4 = scatter(siteXpos(site)+[0:3], MLstats{site,4}(2,:), '^', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth, "MarkerFaceColor", chartColours(site,:), "MarkerFaceAlpha", 0.3);
rec5_3.DisplayName = 'Site 5 Recall';

% Site 6
site = 6;
pre6_3 = scatter(siteXpos(site)+[0:3], MLstats{site,3}(1,:), 'o', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth);
pre6_2 = scatter(siteXpos(site)+[0:3], MLstats{site,2}(1,:), 'v', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth);
pre6_4 = scatter(siteXpos(site)+[0:3], MLstats{site,4}(1,:), '^', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth);
pre6_3.DisplayName = 'Site 6 Precision';
pre6_2.DisplayName = 'Site 6 (0.25)';
pre6_4.DisplayName = 'Site 6 (0.75)';
rec6_3 = scatter(siteXpos(site)+[0:3], MLstats{site,3}(2,:), 'o', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth, "MarkerFaceColor", chartColours(site,:), "MarkerFaceAlpha", 0.3);
rec6_2 = scatter(siteXpos(site)+[0:3], MLstats{site,2}(2,:), 'v', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth, "MarkerFaceColor", chartColours(site,:), "MarkerFaceAlpha", 0.3);
rec6_4 = scatter(siteXpos(site)+[0:3], MLstats{site,4}(2,:), '^', 'MarkerEdgeColor', chartColours(site,:), 'SizeData', 125, "MarkerEdgeAlpha", 1, 'LineWidth', LWidth, "MarkerFaceColor", chartColours(site,:), "MarkerFaceAlpha", 0.3);
rec6_3.DisplayName = 'Site 6 Recall';


xlim(Xlims)
xticks([1:4])
xticklabels({'Field','Vert','45°','Horiz'})
ylim([-0.05, 1.07])
ax1.XGrid = 'on';
ax1.YGrid = 'on';
ax1.FontSize = 11;
title('BirdNET - Precision and Recall', 'FontSize', 13)
set(ax1,'Box','on')


topLG = legend([pre1_3 pre2_3 pre3_3 pre4_3 pre5_3 pre6_3 rec1_3 rec2_3 rec3_3 rec4_3 rec5_3 rec6_3], ...
    'Orientation','Horizontal','NumColumns',6,'FontSize',11);
topLG.Location = 'northoutside';

ax2 = nexttile(2); % BN number of species ---------------------------------
topBar = bar(nSpecsMultThresh{3,1}', 1);
xlim(Xlims)
xticks([1:4])
xticklabels({'Field','Vert','45°','Horiz'})
ax2.FontSize = 11;
ylim([0 6.5])
yticks([0:2:6])
title('BirdNET - No. Species', 'FontSize', 13)
hold on

negErrBarT = nSpecsMultThresh{4,1} - nSpecsMultThresh{3,1};
posErrBarT = nSpecsMultThresh{2,1} - nSpecsMultThresh{3,1};

for k = 1:6 % Display errorbars and number for each bar above bar
    xtips = topBar(k).XEndPoints;
    eb = errorbar(xtips, nSpecsMultThresh{3,1}(k,:), negErrBarT(k,:),...
        posErrBarT(k,:), 'Color', 'black', 'LineWidth', 1);
    eb.LineStyle = 'none';
    labelsT = string(topBar(k).YData);

    text(xtips(posErrBarT(k,:) > 0)-0.03,topBar(k).YEndPoints(...
        posErrBarT(k,:) > 0),labelsT(posErrBarT(k,:) > 0),...
        'HorizontalAlignment','center','VerticalAlignment','bottom')

    text(xtips(posErrBarT(k,:) == 0),nSpecsMultThresh{2,1}(k,...
        posErrBarT(k,:) == 0),labelsT(posErrBarT(k,:) == 0),...
        'HorizontalAlignment','center','VerticalAlignment','bottom')

    topBar(k).FaceColor = chartColours(k, :);
    topBar(k).DisplayName = strcat("Site",{' '},num2str(k));
end

ax2.YGrid = 'on';

barLG = legend([topBar(1) topBar(2) topBar(3) topBar(4) topBar(5) topBar(6)], ...
    'NumColumns',2,'FontSize',11);
barLG.Location = 'northwest';


ax3 = nexttile(3);  % BN number of calls (3 s windows) --------------------
midBar = bar(nCallsMultThresh{3,1}', 1);
xlim(Xlims)
xticks([1:4])
xticklabels({'Field','Vert','45°','Horiz'})
ax3.FontSize = 11;
ylim([0 35])
yticks([0:10:40])
title('BirdNET - No. Predictions', 'FontSize', 13)
hold on

negErrBarM = nCallsMultThresh{4,1} - nCallsMultThresh{3,1};
posErrBarM = nCallsMultThresh{2,1} - nCallsMultThresh{3,1};

for k = 1:6 % Display errorbars and number for each bar above bar
    xtips = midBar(k).XEndPoints;
    eb = errorbar(xtips, nCallsMultThresh{3,1}(k,:), negErrBarM(k,:),...
    posErrBarM(k,:), 'Color', 'black', 'LineWidth', 1);
    eb.LineStyle = 'none';
    ytips = midBar(k).YEndPoints;
    labelsM = string(midBar(k).YData);

    text(xtips(posErrBarM(k,:) > 0 & ytips < 10)-0.03,ytips(...
        posErrBarM(k,:) > 0 & ytips < 10),labelsM(posErrBarM(k,:)...
        > 0 & ytips < 10),...
        'HorizontalAlignment','center','VerticalAlignment','bottom')

    text(xtips(posErrBarM(k,:) > 0 & ytips > 10)-0.05,ytips(...
        posErrBarM(k,:) > 0 & ytips > 10),labelsM(posErrBarM(k,:)...
        > 0 & ytips > 10),...
        'HorizontalAlignment','center','VerticalAlignment','bottom')

    text(xtips(posErrBarM(k,:) == 0),nCallsMultThresh{2,1}(k,...
        posErrBarM(k,:) == 0),labelsM(posErrBarM(k,:) == 0),...
        'HorizontalAlignment','center','VerticalAlignment','bottom')

    midBar(k).FaceColor = chartColours(k, :);
end

ax3.YGrid = 'on';

ax4 = nexttile(4);  % Ornothilogist number of species and calls -----------
lowBar = bar([numManSpecs numManDets], 1);
xticks([1:6])
xticklabels({'Site 1','Site 2','Site 3','Site 4','Site 5','Site 6'})
ax4.FontSize = 11;
ylim([0 170])
yticks([0:50:150])
title('Ornithologist - No. Species and Predictions', 'FontSize', 13)

for k = 1:2 % Display number for each bar above bar
    xtips = lowBar(k).XEndPoints;
    ytips = lowBar(k).YEndPoints;
    labelsB = string(lowBar(k).YData);
    text(xtips,ytips,labelsB,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
    lowBar(k).FaceColor = chartColours(k+6, :);
end

lowBar(1).DisplayName = 'No. Species';
lowBar(2).DisplayName = 'No. Predictions';

lowLG = legend([lowBar(1) lowBar(2)], 'FontSize',11);

