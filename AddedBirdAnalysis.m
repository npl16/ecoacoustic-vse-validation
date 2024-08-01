%% About
% Script to analyse performance of BirdNET and HARKBird at classifying and 
% localising bird sounds added to VSE-replicated soundscapes by two main 
% methods (to be compared): encoded in the Ambisonics domain, and 
% reproduced from directly from individual loudspeakers.

% V3.1, 02.05.2024

% This script is sctructured as follows:

% Section 1. Set initial variables (incl. ground truth data)
% Section 2. Import and process HARKBird data to obtain median
%            localisation/errors for each added bird call
% Section 3. Import and process BirdNET data to obtain mean classification
%            confidence for each added bird call
% Section 4. Calculate overall means and standard deviations of BirdNET 
%            confidence values and HARKBird errors.
% Section 5. Make Figure 7 - 2-pane plot of BirdNET confidence against 
%            added birds' elevation and HARKBird error against elevation
% Section 6. Supp. Figure 2 - plot of HARKBird error against confidence for
%            the 4 methods of adding bird calls explored (ambisonic encoding 
%            with low pass, individual loudspeaker with low pass,
%            individual loudspeaker with no reverb, individual loudspeaker
%            with no soundscape or reverb). 
% Sections 7 & 8. Mann-Whiteney U Tests to assess for significant
%            differences between BirdNET confidence and HARKBird error 
%            between bird calls added by ambisonic encoding and individual 
%            loudspeaker playback.

% Abbreviations:
% ambi = ambisonic(s);
% Az = Azimuth;
% BN = BirdNET;
% conf = BirdNET confidence;
% El = Elevation;
% HB = HARKBird; 
% locs = Localisations (estimated azimuth) of avian calls in HARKBird; 
% LP = low passed; 
% lspk = loudspeaker(s); 
% NR = no reverb; 
% solo = Soloed playback of a bird callfrom an individual loudspeaker 
%        (no background soundscape orreverberation)


%--------------------------------------------------------------------------
%% 1. Create main variables and stuctures (incl. ground truth data) need for importing/analysing data

numABRecs = 10; % Number of recordings with added bird sounds

% Common names of added birds in the order they appear in the recordings
% (ordered first by site and then by the sweep used to reproduce
% reverberation).
addedBirdNames = {'Robin';'Sparrow';'BlueTit';'Pigeon';'Greenfinch';...
    'BlackBird';'Greenfinch';'BlueTit';'Greenfinch';'Robin'};

% Site number and sweep used for reproducing the reverberation of each
% added bird
addedBirdSiteSweep = {'s1s2';'s1s3';'s2s1';'s2s3';'s3s2';'s3s4';'s4s1';'s4s3';'s5s1';'s5s4'};

Site = [1; 1; 2; 2; 3; 3; 4; 4; 5; 5];
Sweep = [2; 3; 1; 3; 2; 4; 1; 3; 1; 4];

ambiAzimuths = [102; 182; 30; 180; 63; 260; 0; 150; 13; 243]; % Azimuths of birds added by ambisonic encoding
ambiElevations = [56; 67; 45; 45; 35; 71; 0; 25; 68; 83]; % Elevations of birds added by ambisonic encoding
adjAmbiAzimuths = ambiAzimuths;
adjAmbiAzimuths(adjAmbiAzimuths > 180) = adjAmbiAzimuths(adjAmbiAzimuths > 180) - 360; % Azimuths converted to -180° to 180° scale

lspkAzimuths = [126; 198; 18; 198; 90; 270; 0; 162; 342; 270]; % Azimuths of birds added by playback from individual loudspeakers
lspkElevations = [56; 56; 32; 56; 32; 56; 0; 32; 56; 56]; % Elevations of birds added by playback from individual loudspeakers
lspks = [28; 29; 21; 29; 22; 30; 11; 23; 26; 30]; % Number of the loudspeaker used (based on internal numbering of louspeakers in VSE)
adjLspkAzimuths = lspkAzimuths;
adjLspkAzimuths(adjLspkAzimuths > 180) = adjLspkAzimuths(adjLspkAzimuths > 180) - 360; % Azimuths converted to -180° to 180° scale

startTimes = [300; 300; 270; 510; 105; 105; 300; 180; 90; 90]; % Start and end times of added bird calls
endTimes = [303; 453; 365; 540; 180; 240; 450; 277; 166; 93];
soloStartTimes = zeros(numABRecs,1); % 'Soloed' birds (with no reverb or background soundscape) start at 0s...
soloEndTimes = endTimes - startTimes; % ... and end after the duration of the added bird call recording.

realAzimuths = [ambiAzimuths, lspkAzimuths, lspkAzimuths, lspkAzimuths];
realElevations = [ambiElevations, lspkElevations, lspkElevations, lspkElevations];

% Create a table with the above details 
addedBirdDetails = table(addedBirdSiteSweep,addedBirdNames,startTimes,...
    endTimes,ambiAzimuths,ambiElevations,lspkAzimuths,lspkElevations,...
    lspks,Site,Sweep,'VariableNames',["SiteAndSweep","CommonName",...
    "StartTime","EndTime","AmbiAzimuth","AmbiElevations","LspkAzimuths",...
    "LspkElevations","LspkNumber","Site","Sweep"]);

addedBirdDetails = convertvars(addedBirdDetails,["SiteAndSweep","CommonName"],"categorical");


%--------------------------------------------------------------------------
%% 2A. Import the HARKBird Data

HBoffset = 30; % Offset between HB's 0 azimuth and the 6-mic's 0 azimuth


allLocs = cell(1,4); % Cell to contain tables of relevant localisations from each of the 4 embedding techniques

for i = 1:4
    for j = 1:numABRecs

        %Import the HB results CSV files
        if i == 1 % For calls added by ambisonic encoding and low-passed
            currentLocs = readtable(strcat(addedBirdSiteSweep{j},addedBirdNames{j},'AmbiLP-HBdets.csv'));
        elseif i == 2 % For calls added by individual loudspeakers and low-passed
            currentLocs = readtable(strcat(addedBirdSiteSweep{j},addedBirdNames{j},'LSPKLP-HBdets.csv'));
        elseif i == 3 % For calls added by individual loudspeakers with no reverb
            currentLocs = readtable(strcat(addedBirdSiteSweep{j},addedBirdNames{j},'LSPKNR-HBdets.csv'));
        else % For calls added by individual loudspeakers with no reverb or soundscape
            currentLocs = readtable(strcat(addedBirdSiteSweep{j},addedBirdNames{j},'Solo-HBdets.csv'));            
        end

        currentLocs.StartTime = floor(currentLocs.StartTime); % Round detections' start times DOWN to nearest second
        currentLocs.StartTime(currentLocs.StartTime<0) = 0; % Make any negative start times 0 (almost never happens, but is an occasional glitch in the HB outputs).
        currentLocs.EndTime = ceil(currentLocs.EndTime); % Round detections' end times UP to nearest second
    
        % Remove HB offset (error between where it thinks 0° is, and the real
        % position) from start and end azimuth values
        currentLocs.StartAzimuth = currentLocs.StartAzimuth - HBoffset; 
        currentLocs.EndAzimuth = currentLocs.EndAzimuth - HBoffset;
    
        % Shift from -180° to 180° scale to 0° to 360° by adding 360° to all
        % negative azimuth values
        currentLocs.StartAzimuth(currentLocs.StartAzimuth<0) = currentLocs.StartAzimuth(currentLocs.StartAzimuth<0) + 360;
        currentLocs.EndAzimuth(currentLocs.EndAzimuth<0) = currentLocs.EndAzimuth(currentLocs.EndAzimuth<0) + 360;
    
        minTime = min(currentLocs.StartTime);
        maxTime = max(currentLocs.EndTime);
    
        numLocsPerSec = zeros(maxTime,1); % Store the number of locs in each second
    
        for p = 1:height(currentLocs) % For every entry in the current locs
            for q = (currentLocs.StartTime(p)+1):currentLocs.EndTime(p) % For every second in each entry
                numLocsPerSec(q,1) = numLocsPerSec(q,1) + 1; % tally the number of times those seconds appear
            end
        end
        
        % Create matrix to store localisations for each second (there are 
        % often multiple per second); initially all values set to 370 to 
        % identify seconds with fewer than this recording's max. number of
        % locs in 1 second:
        locsPerSec = 370*ones(maxTime,max(numLocsPerSec)); 

        for p = 1:height(currentLocs) % For every entry in the current locs
            for q = (currentLocs.StartTime(p)+1):currentLocs.EndTime(p) % For every second in each entry
                
                k = 1; % Initialise a counter k
                    
                    while k <= numLocsPerSec(q,1) % While k is less than the number of times this second appears in the current HB results
                        if locsPerSec(q,k) ~= 370 % If the Locations column for this second is not 370 at the current k (i.e. if the 'slot' in the current column is already taken by a previous detection of this column)...
                            k = k + 1; % Increase k (until we hit a column for this second that = 370 (i.e. a column for this second that is free))
                        
                        else
                            break
                        end
                    end
                    
                    % Once we have a free 'slot' for this second, add the average
                    % azimuth value here:
                    locsPerSec(q,k) = mean([currentLocs.StartAzimuth(p),currentLocs.EndAzimuth(p)]);        
            end
        end
    
        % Filter to relevant start and end times for added bird:
        if i < 4
            locsPerSec = locsPerSec(startTimes(j)+1:endTimes(j), :);
        else
            if soloEndTimes(j) < maxTime
                locsPerSec = locsPerSec(soloStartTimes(j)+1:soloEndTimes(j), :);
            end
        end
    
        currentAvLocs = 370*ones(height(locsPerSec),1); % Store array of averaged Locs for each second for current rec here
    
        for a = 1:height(locsPerSec) % For each second
    
            if locsPerSec(a,1) ~= 370 % If there are locs at this second
    
                % First get the magnitude of the azimuths' shortest
                % angle from 0°
                LocsToCheck = locsPerSec(a,find(locsPerSec(a,:) < 370)); % Only consider non-370 values (i.e. actual HB locs)
                
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
                currentAvLocs(currentAvLocs(a,1) < 0, 1) = currentAvLocs(currentAvLocs(a,1) < 0, 1) + 360; % Convert any negative values back to a 0-360° scale
            end
        end

        currentAvLocs(currentAvLocs < 0, 1) = currentAvLocs(currentAvLocs < 0, 1) + 360; % Convert any negative values back to a 0-360° scale


        % Put the data into tables:

        numSecs = height(currentAvLocs);
        
        if i < 4
            % For the first 3 conditions
            currentHBres = table(Site(j)*ones(numSecs,1), Sweep(j)*ones(numSecs,1),...
                [startTimes(j):endTimes(j)-1]', repmat(addedBirdNames(j), numSecs, 1),...
                currentAvLocs, realAzimuths(j,i)*ones(numSecs,1), ...
                realElevations(j,i)*ones(numSecs,1), 'VariableNames',...
                ["Site","Sweep","Time","CommonName","HB_Azimuth","RealAzimuth","RealElevation"]);

        else 
            % For soloed bird calls (as these have different time stamps):
            currentHBres = table(Site(j)*ones(numSecs,1), Sweep(j)*ones(numSecs,1),...
                [soloStartTimes(j):numSecs-1]', repmat(addedBirdNames(j), numSecs, 1),...
                currentAvLocs, realAzimuths(j,i)*ones(numSecs,1), ...
                realElevations(j,i)*ones(numSecs,1), 'VariableNames',...
                ["Site","Sweep","Time","CommonName","HB_Azimuth","RealAzimuth","RealElevation"]);
        end

        allLocs{1,i} = [allLocs{1,i}; currentHBres]; % add to allLocs table
        allLocs{1,i} = convertvars(allLocs{1,i},["CommonName"],"categorical");
        

    end
end

%Extract the tables of all 4 embedding methods' locs from allLocs and add
%to new tables with additional columns for error etc.
ambiLocs = allLocs{1,1}(allLocs{1,1}.HB_Azimuth < 370,:); % Only include actual results (i.e. with azimuth values < 370)
ambiLocs.Error =  ambiLocs.HB_Azimuth - ambiLocs.RealAzimuth; % Calculate error 
ambiLocs.AdjError = ambiLocs.Error;
ambiLocs.AdjError(ambiLocs.AdjError <= -180) = ambiLocs.AdjError(ambiLocs.AdjError <= -180) + 360; % Create additional column where error is expressed on a 0-360 scale
ambiLocs.AdjError(ambiLocs.AdjError > 180) = ambiLocs.AdjError(ambiLocs.AdjError > 180) - 360;
ambiLocs.AdjRealAzimuth = ambiLocs.RealAzimuth;
ambiLocs.AdjRealAzimuth(ambiLocs.AdjRealAzimuth > 180) = ambiLocs.AdjRealAzimuth(ambiLocs.AdjRealAzimuth > 180) - 360; % Create additional column where the actual azimuth of added birds is expressed on -180 to 180 scale

%Repeat for other encoding methods
lspkLPLocs = allLocs{1,2}(allLocs{1,2}.HB_Azimuth < 370,:);
lspkLPLocs.Error = lspkLPLocs.HB_Azimuth - lspkLPLocs.RealAzimuth;
lspkLPLocs.AdjError = lspkLPLocs.Error;
lspkLPLocs.AdjError(lspkLPLocs.AdjError <= -180) = lspkLPLocs.AdjError(lspkLPLocs.AdjError <= -180) + 360;
lspkLPLocs.AdjError(lspkLPLocs.AdjError > 180) = lspkLPLocs.AdjError(lspkLPLocs.AdjError > 180) - 360;
lspkLPLocs.AdjRealAzimuth = lspkLPLocs.RealAzimuth;
lspkLPLocs.AdjRealAzimuth(lspkLPLocs.AdjRealAzimuth > 180) = lspkLPLocs.AdjRealAzimuth(lspkLPLocs.AdjRealAzimuth > 180) - 360;

lspkNRLocs = allLocs{1,3}(allLocs{1,3}.HB_Azimuth < 370,:);
lspkNRLocs.Error = lspkNRLocs.HB_Azimuth - lspkNRLocs.RealAzimuth;
lspkNRLocs.AdjError = lspkNRLocs.Error;
lspkNRLocs.AdjError(lspkNRLocs.AdjError <= -180) = lspkNRLocs.AdjError(lspkNRLocs.AdjError <= -180) + 360;
lspkNRLocs.AdjError(lspkNRLocs.AdjError > 180) = lspkNRLocs.AdjError(lspkNRLocs.AdjError > 180) - 360;
lspkNRLocs.AdjRealAzimuth = lspkNRLocs.RealAzimuth;
lspkNRLocs.AdjRealAzimuth(lspkNRLocs.AdjRealAzimuth > 180) = lspkNRLocs.AdjRealAzimuth(lspkNRLocs.AdjRealAzimuth > 180) - 360;

lspkSoloLocs = allLocs{1,4}(allLocs{1,4}.HB_Azimuth < 370,:);
lspkSoloLocs.Error = lspkSoloLocs.HB_Azimuth - lspkSoloLocs.RealAzimuth;
lspkSoloLocs.AdjError = lspkSoloLocs.Error;
lspkSoloLocs.AdjError(lspkSoloLocs.AdjError <= -180) = lspkSoloLocs.AdjError(lspkSoloLocs.AdjError <= -180) + 360;
lspkSoloLocs.AdjError(lspkSoloLocs.AdjError > 180) = lspkSoloLocs.AdjError(lspkSoloLocs.AdjError > 180) - 360;
lspkSoloLocs.AdjRealAzimuth = lspkSoloLocs.RealAzimuth;
lspkSoloLocs.AdjRealAzimuth(lspkSoloLocs.AdjRealAzimuth > 180) = lspkSoloLocs.AdjRealAzimuth(lspkSoloLocs.AdjRealAzimuth > 180) - 360;


%--------------------------------------------------------------------------
%% 2B. Create tables with summary statistics (means and medians of HARKBird errors)

% For each embedding method, create new empty table to contain 1 row for 
% each added bird alongside empty arrays to store the means and medians of 
% HARKBird localisation errors for each added bird:
avAmbiLocs = table;
avAmbiError = 370*ones(10,1);
avAmbiAdjError = 370*ones(10,1);
medAmbiError = 370*ones(10,1);
medAmbiAdjError = 370*ones(10,1);

avLspkLPLocs = table;
avLspkLPError = 370*ones(10,1);
avLspkLPAdjError = 370*ones(10,1);
medLspkLPError = 370*ones(10,1);
medLspkLPAdjError = 370*ones(10,1);

avLspkNRLocs = table;
avLspkNRError = 370*ones(10,1);
avLspkNRAdjError = 370*ones(10,1);
medLspkNRError = 370*ones(10,1);
medLspkNRAdjError = 370*ones(10,1);

avLspkSoloLocs = table;
avLspkSoloError = 370*ones(10,1);
avLspkSoloAdjError = 370*ones(10,1);
medLspkSoloError = 370*ones(10,1);
medLspkSoloAdjError = 370*ones(10,1);

% Create empty arrays to store HB errors with outliers removed
ambiLocAdjErrors = []; % Columns will be: Site, Sweep, AdjError
lspkLPLocAdjErrors = [];
lspkNRLocAdjErrors = [];
lspkSoloLocAdjErrors = [];


for a = 1:10 % For each added bird

    % For ambisonically-encoded bird sounds
    % Select results from previous table for each individual bird
    newLocsAmbi = ambiLocs(ambiLocs.Site == Site(a),:);
    newLocsAmbi = newLocsAmbi(newLocsAmbi.Sweep == Sweep(a),:);

    if height(newLocsAmbi) == 0
        % If the bird was not detected, add it to the summary table with 
        % dummy value of '370' for its azimuth error
        newLocsAmbi = table(Site(a),Sweep(a),startTimes(a),addedBirdNames(a),...
            370,ambiAzimuths(a),ambiElevations(a),'VariableNames',...
            ["Site","Sweep","Time","CommonName","HB_Azimuth","RealAzimuth", "RealElevation"]);
    
    else
        % Calculate mean and median of HARKBird error (with outliers
        % removed) for ambisonically-encoded low-passed calls
        avAmbiError(a) = mean(rmoutliers(newLocsAmbi.Error));
        noOutlierAmbiAdjError = rmoutliers(newLocsAmbi.AdjError);
        avAmbiAdjError(a) = mean(noOutlierAmbiAdjError);
        
        ambiLocAdjErrors = [ambiLocAdjErrors; Site(a)*...
            ones(height(noOutlierAmbiAdjError),1), ...
            Sweep(a)*ones(height(noOutlierAmbiAdjError),1), ...
            noOutlierAmbiAdjError];

        medAmbiError(a) = median(rmoutliers(newLocsAmbi.Error));
        medAmbiAdjError(a) = median(noOutlierAmbiAdjError);
    end

    avAmbiLocs = [avAmbiLocs; newLocsAmbi(1,[1:7])];


    % Repeat for direct-from-loudspeaker low-passed calls
    newLocsLspk = lspkLPLocs(lspkLPLocs.Site == Site(a),:);
    newLocsLspk = newLocsLspk(newLocsLspk.Sweep == Sweep(a),:);

    if height(newLocsLspk) == 0
        newLocsLspk = table(Site(a),Sweep(a),startTimes(a),addedBirdNames(a),...
            370,lspkAzimuths(a),lspkElevations(a),'VariableNames',...
            ["Site","Sweep","Time","CommonName","HB_Azimuth","RealAzimuth", "RealElevation"]);
    else
        avLspkLPError(a) = mean(rmoutliers(newLocsLspk.Error));
        noOutlierLspkAdjError = rmoutliers(newLocsLspk.AdjError);
        avLspkLPAdjError(a) = mean(noOutlierLspkAdjError);

        lspkLPLocAdjErrors = [lspkLPLocAdjErrors; Site(a)*...
            ones(height(noOutlierLspkAdjError),1), ...
            Sweep(a)*ones(height(noOutlierLspkAdjError),1), ...
            noOutlierLspkAdjError];

        medLspkLPError(a) = median(rmoutliers(newLocsLspk.Error));
        medLspkLPAdjError(a) = median(rmoutliers(newLocsLspk.AdjError));
    end

    avLspkLPLocs = [avLspkLPLocs; newLocsLspk(1,[1:7])];


    % For direct-from-loudspeaker no-reverb calls
    newLocsLspk = lspkNRLocs(lspkNRLocs.Site == Site(a),:);
    newLocsLspk = newLocsLspk(newLocsLspk.Sweep == Sweep(a),:);

    if height(newLocsLspk) == 0
        newLocsLspk = table(Site(a),Sweep(a),startTimes(a),addedBirdNames(a),...
            370,lspkAzimuths(a),lspkElevations(a),'VariableNames',...
            ["Site","Sweep","Time","CommonName","HB_Azimuth","RealAzimuth", "RealElevation"]);
    else
        avLspkNRError(a) = mean(rmoutliers(newLocsLspk.Error));
        noOutlierLspkAdjError = rmoutliers(newLocsLspk.AdjError);
        avLspkNRAdjError(a) = mean(noOutlierLspkAdjError);

        lspkNRLocAdjErrors = [lspkNRLocAdjErrors; Site(a)*...
            ones(height(noOutlierLspkAdjError),1), ...
            Sweep(a)*ones(height(noOutlierLspkAdjError),1), ...
            noOutlierLspkAdjError];

        medLspkNRError(a) = median(rmoutliers(newLocsLspk.Error));
        medLspkNRAdjError(a) = median(rmoutliers(newLocsLspk.AdjError));
    end

    avLspkNRLocs = [avLspkNRLocs; newLocsLspk(1,[1:7])];

    % Repeat for direct-from-loudspeaker no-reverb calls
    newLocsLspk = lspkSoloLocs(lspkSoloLocs.Site == Site(a),:);
    newLocsLspk = newLocsLspk(newLocsLspk.Sweep == Sweep(a),:);

    if height(newLocsLspk) == 0
        newLocsLspk = table(Site(a),Sweep(a),startTimes(a),addedBirdNames(a),...
            370,lspkAzimuths(a),lspkElevations(a),'VariableNames',...
            ["Site","Sweep","Time","CommonName","HB_Azimuth","RealAzimuth", "RealElevation"]);
    else
        avLspkSoloError(a) = mean(rmoutliers(newLocsLspk.Error));
        noOutlierLspkAdjError = rmoutliers(newLocsLspk.AdjError);
        avLspkSoloAdjError(a) = mean(noOutlierLspkAdjError);

        lspkSoloLocAdjErrors = [lspkSoloLocAdjErrors; Site(a)*...
            ones(height(noOutlierLspkAdjError),1), ...
            Sweep(a)*ones(height(noOutlierLspkAdjError),1), ...
            noOutlierLspkAdjError];

        medLspkSoloError(a) = median(rmoutliers(newLocsLspk.Error));
        medLspkSoloAdjError(a) = median(rmoutliers(newLocsLspk.AdjError));
    end

    avLspkSoloLocs = [avLspkSoloLocs; newLocsLspk(1,[1:7])];

end

% Summary of HARKBird errors for ambisonically-encoded birds
avAmbiLocs.AvError_Ambi = avAmbiError;
avAmbiLocs.AvAdjError_Ambi = avAmbiAdjError;
avAmbiLocs.MedError_Ambi = medAmbiError;
avAmbiLocs.MedAdjError_Ambi = medAmbiAdjError;

% Summary of HARKBird errors for direct-from-loudspeaker low-passed birds
avLspkLPLocs.AvError_LspkLP = avLspkLPError;
avLspkLPLocs.AvAdjError_LspkLP = avLspkLPAdjError;
avLspkLPLocs.MedError_LspkLP = medLspkLPError;
avLspkLPLocs.MedAdjError_LspkLP = medLspkLPAdjError;

% Summary of HARKBird errors for direct-from-loudspeaker no-reverb birds
avLspkNRLocs.AvError_LspkNR = avLspkNRError;
avLspkNRLocs.AvAdjError_LspkNR = avLspkNRAdjError;
avLspkNRLocs.MedError_LspkNR = medLspkNRError;
avLspkNRLocs.MedAdjError_LspkNR = medLspkNRAdjError;

% Summary of HARKBird errors for direct-from-loudspeaker soloed birds
avLspkSoloLocs.AvError_LspkSolo = avLspkSoloError;
avLspkSoloLocs.AvAdjError_LspkSolo = avLspkSoloAdjError;
avLspkSoloLocs.MedError_LspkSolo = medLspkSoloError;
avLspkSoloLocs.MedAdjError_LspkSolo = medLspkSoloAdjError;

% Compile all separate tables into one:
avLocErrors = [avAmbiLocs(:,[1,2,3,4,6,7,8,9,10,11]),avLspkLPLocs(:,[8,9,10,11]),avLspkNRLocs(:,[8,9,10,11]),avLspkSoloLocs(:,[8,9,10,11])];

% Create separate table with just the median (rather than mean) errors:
medLocErrors = avLocErrors(:,[1,2,3,4,5,6,10,14,18,22]);



%--------------------------------------------------------------------------
%% 3A. Import the BirdNET data

% Import the overall data

addedBirdsBNambi = readtable('AllAddedBirds-BNResults-Ambi.csv');
addedBirdsBNambi = convertvars(addedBirdsBNambi,["ScientificName","CommonName"],"categorical");

addedBirdsBNlspk = readtable('AllAddedBirds-BNResults-Lspk.csv');
addedBirdsBNlspk = convertvars(addedBirdsBNlspk,["ScientificName","CommonName"],"categorical");

%% 3B. Filter BN data down to each added bird

% Create empty tables and arrays to store mean and median BN confidence for
% each added bird:
avAmbiBNres = table;
avAmbiConfs = 370*ones(10,1);
medAmbiConfs = 370*ones(10,1);

avLspkBNres = table;
avLspkConfs = 370*ones(10,1);
medLspkConfs = 370*ones(10,1);

% Calculate the mean and median BN confidence value for each added bird
for a = 1:10 % For each added bird...
    
    % For ambisonically-encoded bird sounds
    % Store the current added birds BN results in new table
    newConfsAmbi = addedBirdsBNambi(addedBirdsBNambi.Site == Site(a),[1,4,5,6,7]);
    newConfsAmbi = newConfsAmbi(newConfsAmbi.Sweep == Sweep(a),:);

    if height(newConfsAmbi) == 0 % If bird not detected, create dummy entry
        newConfsAmbi = table(startTimes(a),addedBirdNames(a),0,Site(a),Sweep(a),'VariableNames',...
            ["Start_s_","CommonName","Confidence","Site","Sweep"]);
    else % Get mean and median confidence
        avAmbiConfs(a) = mean(newConfsAmbi.Confidence);
        medAmbiConfs(a) = median(newConfsAmbi.Confidence);

    end
    
    % Start to construct summary table - add select details (start 
    % time, bird name, site, sweep) from first row of current table to
    % existing entries:
    avAmbiBNres = [avAmbiBNres; newConfsAmbi(1,[1,2,4,5])];


    % Repeat for direct-from-loudspeaker sounds
    newConfsLspk = addedBirdsBNlspk(addedBirdsBNlspk.Site == Site(a),[1,4,5,6,7]);
    newConfsLspk = newConfsLspk(newConfsLspk.Sweep == Sweep(a),:);

    if height(newConfsLspk) == 0
        newConfsLspk = table(startTimes(a),addedBirdNames(a),0,Site(a),Sweep(a),'VariableNames',...
            ["Start_s_","CommonName","Confidence","Site","Sweep"]);
    else
        avLspkConfs(a) = mean(newConfsLspk.Confidence);
        medLspkConfs(a) = median(newConfsLspk.Confidence);

    end

    avLspkBNres = [avLspkBNres; newConfsLspk(1,[1,2,4,5])];

end

% Add final averaged confidence values to summary tables
avAmbiBNres.AvConfAmbi = avAmbiConfs;
avAmbiBNres.MedConfAmbi = medAmbiConfs;

avLspkBNres.AvConfLspk = avLspkConfs;
avLspkBNres.MedConfLspk = medLspkConfs;

% Combine data into one overall summary table with HB localisaton errors
addedBirdSummary = [medLocErrors,avAmbiBNres(:,5),avLspkBNres(:,5)];


%--------------------------------------------------------------------------
%% 4. Overall means and standard deviations of BN confidence values and HB errors

overallMeanAmbiConfs = mean(addedBirdSummary.AvConfAmbi(addedBirdSummary.AvConfAmbi<370));
overallStDAmbiConfs = std(addedBirdSummary.AvConfAmbi(addedBirdSummary.AvConfAmbi<370));

overallMeanLspkLPConfs = mean(addedBirdSummary.AvConfLspk(addedBirdSummary.AvConfLspk<370));
overallStDLspkLPConfs = std(addedBirdSummary.AvConfLspk(addedBirdSummary.AvConfLspk<370));

overallMeanAmbiLocErrors = mean(abs(addedBirdSummary.MedAdjError_Ambi(addedBirdSummary.MedAdjError_Ambi<370)));
overallStDAmbiLocErrors = std(abs(addedBirdSummary.MedAdjError_Ambi(addedBirdSummary.MedAdjError_Ambi<370)));

overallMeanLspkLPLocErrors = mean(abs(addedBirdSummary.MedAdjError_LspkLP(addedBirdSummary.MedAdjError_LspkLP<370)));
overallStDLspkLPLocErrors = std(abs(addedBirdSummary.MedAdjError_LspkLP(addedBirdSummary.MedAdjError_LspkLP<370)));

overallMeanLspkNRLocErrors = mean(abs(addedBirdSummary.MedAdjError_LspkNR(addedBirdSummary.MedAdjError_LspkNR<370)));
overallStDLspkNRLocErrors = std(abs(addedBirdSummary.MedAdjError_LspkNR(addedBirdSummary.MedAdjError_LspkNR<370)));

overallMeanLspkSoloLocErrors = mean(abs(addedBirdSummary.MedAdjError_LspkSolo(addedBirdSummary.MedAdjError_LspkSolo<370)));
overallStDLspkSoloLocErrors = std(abs(addedBirdSummary.MedAdjError_LspkSolo(addedBirdSummary.MedAdjError_LspkSolo<370)));


%--------------------------------------------------------------------------
%% 5. Fig. 7 - Plot of added bird true elevation against BN confidence and HB error 
% 2-by-1 tiled plot with median absolute HB error against true elevation
% and mean BN confidence against true elevation.
% Only Ambi and Direct from Loudspeaker (LP) encoding

figure
tiledlayout(1,2,'TileSpacing','loose','Padding','compact') % Create 2-tile figure

currentFig = gcf;
set(currentFig, 'Units', 'centimeters', 'Position', [0, 0, 25, 22], 'PaperUnits', 'centimeters', 'PaperSize', [25, 22]) % Set figure dimensions


% Plot of median absolute HB error against true elevation with different 
% markers per species:

ax1 = nexttile(1);
% Plot each species with a different marker (for mapping, see Display Names
% for objects c1 to c12 below)
h1 = scatter(abs(ambiElevations([1,10])),abs(medLocErrors.MedAdjError_Ambi([1,10])),100,'vb','LineWidth',1.5);
hold on
grid on
h2 = scatter(abs(ambiElevations(2)),abs(medLocErrors.MedAdjError_Ambi(2)),100,'ob','LineWidth',1.5);
h3 = scatter(abs(ambiElevations([3,8])),abs(medLocErrors.MedAdjError_Ambi([3,8])),100,'^b','LineWidth',1.5);
h4 = scatter(abs(ambiElevations(4)),abs(medLocErrors.MedAdjError_Ambi(4)),100,'+b','LineWidth',1.5);
h5 = scatter(abs(ambiElevations([5,9])),abs(medLocErrors.MedAdjError_Ambi([5,9])),100,'xb','LineWidth',1.5);
h6 = scatter(abs(ambiElevations([6,7])),abs(medLocErrors.MedAdjError_Ambi([6,7])),100,'squareb','LineWidth',1.5);

h7 = scatter(abs(lspkElevations([1,10])),abs(medLocErrors.MedAdjError_LspkLP([1,10])),100,'vr','LineWidth',1.5);
h8 = scatter(abs(lspkElevations(2)),abs(medLocErrors.MedAdjError_LspkLP(2)),100,'or','LineWidth',1.5);
h9 = scatter(abs(lspkElevations([3,8])),abs(medLocErrors.MedAdjError_LspkLP([3,8])),100,'^r','LineWidth',1.5);
h10 = scatter(abs(lspkElevations(4)),abs(medLocErrors.MedAdjError_LspkLP(4)),100,'+r','LineWidth',1.5);
h11 = scatter(abs(lspkElevations([5,9])),abs(medLocErrors.MedAdjError_LspkLP([5,9])),100,'xr','LineWidth',1.5);
h12 = scatter(abs(lspkElevations([6,7])),abs(medLocErrors.MedAdjError_LspkLP([6,7])),100,'squarer','LineWidth',1.5);


axis([0 90 0 180]) % Set axes' limits
xticks([0:10:90]) % Set steps for x-axis
ax1.FontSize = 11;
ylabel('Absolute Error (deg)','FontSize',14)
xlabel('Added Bird True Elevation (deg)','FontSize',14)
title({'Absolute HARKBird Error against', 'True Elevation'},'FontSize',14)

% Use polyfit and polyval to fit linear model to HB error
[p_ambiElAbsErr,S_ambiElAbsErr] = polyfit(ambiElevations,abs(medLocErrors.MedAdjError_Ambi),1);
fAmbi_ElAbs = polyval(p_ambiElAbsErr,[0:90]);
h25 = plot([0:90],fAmbi_ElAbs,'-b','LineWidth',1);
h25.DisplayName = 'Linear Fit';

[p_lspkLPElAbsErr,S_lspkLPElAbsErr] = polyfit(lspkElevations(1:9),abs(medLocErrors.MedAdjError_LspkLP(1:9)),1);
fLspkLP_ElAbs = polyval(p_lspkLPElAbsErr,[0:90]);
h26 = plot([0:90],fLspkLP_ElAbs,'-r','LineWidth',1);
h26.DisplayName = 'Linear Fit';


hold off


% BN Confidence against elevation
ax2 = nexttile(2);
c1 = scatter(abs(ambiElevations([1,10])),abs(addedBirdSummary.AvConfAmbi([1,10])),100,'vb','LineWidth',1.5);
hold on
grid on
c2 = scatter(abs(ambiElevations(2)),abs(addedBirdSummary.AvConfAmbi(2)),100,'ob','LineWidth',1.5);
c3 = scatter(abs(ambiElevations([3,8])),abs(addedBirdSummary.AvConfAmbi([3,8])),100,'^b','LineWidth',1.5);
c4 = scatter(abs(ambiElevations(4)),abs(addedBirdSummary.AvConfAmbi(4)),100,'+b','LineWidth',1.5);
c5 = scatter(abs(ambiElevations([5,9])),abs(addedBirdSummary.AvConfAmbi([5,9])),100,'xb','LineWidth',1.5);
c6 = scatter(abs(ambiElevations([6,7])),abs(addedBirdSummary.AvConfAmbi([6,7])),100,'squareb','LineWidth',1.5);

c7 = scatter(abs(lspkElevations([1,10])),abs(addedBirdSummary.AvConfLspk([1,10])),100,'vr','LineWidth',1.5);
c8 = scatter(abs(lspkElevations(2)),abs(addedBirdSummary.AvConfLspk(2)),100,'or','LineWidth',1.5);
c9 = scatter(abs(lspkElevations([3,8])),abs(addedBirdSummary.AvConfLspk([3,8])),100,'^r','LineWidth',1.5);
c10 = scatter(abs(lspkElevations(4)),abs(addedBirdSummary.AvConfLspk(4)),100,'+r','LineWidth',1.5);
c11 = scatter(abs(lspkElevations([5,9])),abs(addedBirdSummary.AvConfLspk([5,9])),100,'xr','LineWidth',1.5);
c12 = scatter(abs(lspkElevations([6,7])),abs(addedBirdSummary.AvConfLspk([6,7])),100,'squarer','LineWidth',1.5);


axis([0 90 0.3 1])
xticks([0:10:90])
ax2.FontSize = 11;
ylabel('Confidence','FontSize',14)
xlabel('Added Bird True Elevation (deg)','FontSize',14)
title({'BirdNET Confidence against', 'True Elevation'},'FontSize',14)

% Use polyfit and polyval to fit linear model to BN confidence
[p_ambiElConf,S_ambiElConf] = polyfit(abs(ambiElevations),abs(addedBirdSummary.AvConfAmbi),1);
fAmbiElConf = polyval(p_ambiElConf,[0:90]);
c13 = plot([0:90],fAmbiElConf,'-b','LineWidth',1);
c13.DisplayName = 'Linear Fit';

[p_lspkElConf,S_lspkElConf] = polyfit(abs(lspkElevations([1,2,4,5,6,7,9,10])),abs(addedBirdSummary.AvConfLspk([1,2,4,5,6,7,9,10])),1);
fLspkElConf = polyval(p_lspkElConf,[0:90]);
c14 = plot([0:90],fLspkElConf,'-r','LineWidth',1);
c14.DisplayName = 'Linear Fit';


% Create set of white (invisible) points to create 'titles' within legend
c15 = scatter(45,90,100,'xw','LineWidth',1.5);
c16 = scatter(45,90,100,'ow','LineWidth',1.5);

% Assign disyplay names to plotted points so their names appear in legend
c1.DisplayName = 'Robin';  % assign legend string
c2.DisplayName = 'Sparrow';  % assign legend string
c3.DisplayName = 'Blue Tit';  % assign legend string
c4.DisplayName = 'Pigeon';  % assign legend string
c5.DisplayName = 'Greenfinch';  % assign legend string
c6.DisplayName = 'Blackbird';  % assign legend string
c7.DisplayName = 'Robin';  % assign legend string
c8.DisplayName = 'Sparrow';  % assign legend string
c9.DisplayName = 'Blue Tit';  % assign legend string
c10.DisplayName = 'Pigeon';  % assign legend string
c11.DisplayName = 'Greenfinch';  % assign legend string
c12.DisplayName = 'Blackbird';  % assign legend string


c15.DisplayName = 'AMBISONICS';  % assign legend string
c16.DisplayName = 'FROM LOUDSPEAKER';  % assign legend string


% Create and format legend
lg = legend([c15 c16 ...
    c1 c7  ...
    c2 c8 ...
    c3 c9 ...
    c4 c10 ...
    c5 c11 ...
    c6 c12 ...
    c13 c14],'Orientation','Horizontal','NumColumns',2,'FontSize',11);
lg.NumColumns = 2;
lg.Layout.Tile = 'South';

hold off



%--------------------------------------------------------------------------
%% 6. Supporting Material Fig. 2
% Plot of median absolute HB error against true elevation - All Added Bird 
% Encoding Types (ambiLP, lspkLP, lspkNR, lspkSolo)

figure
%tiledlayout(1,2,'TileSpacing','loose','Padding','compact')

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 27, 25], 'PaperUnits', 'centimeters', 'PaperSize', [27, 25]) % Set dimensions

ax3 = gca;

% For ambisonic encoding with low pass
h1 = scatter(abs(ambiElevations([1,10])),abs(medLocErrors.MedAdjError_Ambi([1,10])),100,'vb','LineWidth',1.5);
hold on
grid on
h2 = scatter(abs(ambiElevations(2)),abs(medLocErrors.MedAdjError_Ambi(2)),100,'ob','LineWidth',1.5);
h3 = scatter(abs(ambiElevations([3,8])),abs(medLocErrors.MedAdjError_Ambi([3,8])),100,'^b','LineWidth',1.5);
h4 = scatter(abs(ambiElevations(4)),abs(medLocErrors.MedAdjError_Ambi(4)),100,'+b','LineWidth',1.5);
h5 = scatter(abs(ambiElevations([5,9])),abs(medLocErrors.MedAdjError_Ambi([5,9])),100,'xb','LineWidth',1.5);
h6 = scatter(abs(ambiElevations([6,7])),abs(medLocErrors.MedAdjError_Ambi([6,7])),100,'squareb','LineWidth',1.5);

% For direct-from-loudspeaker with low pass
h7 = scatter(abs(lspkElevations([1,10])),abs(medLocErrors.MedAdjError_LspkLP([1,10])),100,'vr','LineWidth',1.5);
h8 = scatter(abs(lspkElevations(2)),abs(medLocErrors.MedAdjError_LspkLP(2)),100,'or','LineWidth',1.5);
h9 = scatter(abs(lspkElevations([3,8])),abs(medLocErrors.MedAdjError_LspkLP([3,8])),100,'^r','LineWidth',1.5);
h10 = scatter(abs(lspkElevations(4)),abs(medLocErrors.MedAdjError_LspkLP(4)),100,'+r','LineWidth',1.5);
h11 = scatter(abs(lspkElevations([5,9])),abs(medLocErrors.MedAdjError_LspkLP([5,9])),100,'xr','LineWidth',1.5);
h12 = scatter(abs(lspkElevations([6,7])),abs(medLocErrors.MedAdjError_LspkLP([6,7])),100,'squarer','LineWidth',1.5);

% For direct-from-loudspeaker with no reverb
h13 = scatter(abs(lspkElevations([1,10])),abs(medLocErrors.MedAdjError_LspkNR([1,10])),100,'vm','LineWidth',1.5);
h14 = scatter(abs(lspkElevations(2)),abs(medLocErrors.MedAdjError_LspkNR(2)),100,'om','LineWidth',1.5);
h15 = scatter(abs(lspkElevations([3,8])),abs(medLocErrors.MedAdjError_LspkNR([3,8])),100,'^m','LineWidth',1.5);
h16 = scatter(abs(lspkElevations(4)),abs(medLocErrors.MedAdjError_LspkNR(4)),100,'+m','LineWidth',1.5);
h17 = scatter(abs(lspkElevations([5,9])),abs(medLocErrors.MedAdjError_LspkNR([5,9])),100,'xm','LineWidth',1.5);
h18 = scatter(abs(lspkElevations([6,7])),abs(medLocErrors.MedAdjError_LspkNR([6,7])),100,'squarem','LineWidth',1.5);

% For direct-from-loudspeaker with no reverb or soundscape
h19 = scatter(abs(lspkElevations([1,10])),abs(medLocErrors.MedAdjError_LspkSolo([1,10])),100,'vg','LineWidth',1.5);
h20 = scatter(abs(lspkElevations(2)),abs(medLocErrors.MedAdjError_LspkSolo(2)),100,'og','LineWidth',1.5);
h21 = scatter(abs(lspkElevations([3,8])),abs(medLocErrors.MedAdjError_LspkSolo([3,8])),100,'^g','LineWidth',1.5);
h22 = scatter(abs(lspkElevations(4)),abs(medLocErrors.MedAdjError_LspkSolo(4)),100,'+g','LineWidth',1.5);
h23 = scatter(abs(lspkElevations([5,9])),abs(medLocErrors.MedAdjError_LspkSolo([5,9])),100,'xg','LineWidth',1.5);
h24 = scatter(abs(lspkElevations([6,7])),abs(medLocErrors.MedAdjError_LspkSolo([6,7])),100,'squareg','LineWidth',1.5);

% Set and label axes and title
axis([0 90 0 180])
xticks([0:10:90])
ax3.FontSize = 11;
ylabel('Absolute Error (deg)','FontSize',14)
xlabel('Added Bird True Elevation (deg)','FontSize',14)
title({'Absolute HARKBird Error against True Elevation'},'FontSize',14)

% Fit linear models
[p_ambiElAbsErr,S_ambiElAbsErr] = polyfit(ambiElevations,abs(medLocErrors.MedAdjError_Ambi),1);
fAmbi_ElAbs = polyval(p_ambiElAbsErr,[0:90]);
h25 = plot([0:90],fAmbi_ElAbs,'-b','LineWidth',1);
h25.DisplayName = 'Linear Fit';

[p_lspkLPElAbsErr,S_lspkLPElAbsErr] = polyfit(lspkElevations(1:9),abs(medLocErrors.MedAdjError_LspkLP(1:9)),1);
fLspkLP_ElAbs = polyval(p_lspkLPElAbsErr,[0:90]);
h26 = plot([0:90],fLspkLP_ElAbs,'-r','LineWidth',1);
h26.DisplayName = 'Linear Fit';

[p_lspkNRElAbsErr,S_lspkNRElAbsErr] = polyfit(lspkElevations,abs(medLocErrors.MedAdjError_LspkNR),1);
fLspkNR_ElAbs = polyval(p_lspkNRElAbsErr,[0:90]);
h27 = plot([0:90],fLspkNR_ElAbs,'-m','LineWidth',1);
h27.DisplayName = 'Linear Fit';

[p_lspkSoloElAbsErr,S_lspkSoloElAbsErr] = polyfit(lspkElevations(1:9),abs(medLocErrors.MedAdjError_LspkSolo(1:9)),1);
fLspkSolo_ElAbs = polyval(p_lspkSoloElAbsErr,[0:90]);
h28 = plot([0:90],fLspkSolo_ElAbs,'-g','LineWidth',1);
h28.DisplayName = 'Linear Fit';

% Create set of white (invisible) points to create 'titles' within legend
h29 = scatter(45,90,100,'xw','LineWidth',1.5);
h30 = scatter(45,90,100,'ow','LineWidth',1.5);
h31 = scatter(45,90,100,'+w','LineWidth',1.5);
h32 = scatter(45,90,100,'*w','LineWidth',1.5);

% Assign disyplay names to plotted points so their names appear in legend
h1.DisplayName = 'Robin';  % assign legend string
h2.DisplayName = 'Sparrow';  % assign legend string
h3.DisplayName = 'Blue Tit';  % assign legend string
h4.DisplayName = 'Pigeon';  % assign legend string
h5.DisplayName = 'Greenfinch';  % assign legend string
h6.DisplayName = 'Blackbird';  % assign legend string
h7.DisplayName = 'Robin';  % assign legend string
h8.DisplayName = 'Sparrow';  % assign legend string
h9.DisplayName = 'Blue Tit';  % assign legend string
h10.DisplayName = 'Pigeon';  % assign legend string
h11.DisplayName = 'Greenfinch';  % assign legend string
h12.DisplayName = 'Blackbird';  % assign legend string
h13.DisplayName = 'Robin';  % assign legend string
h14.DisplayName = 'Sparrow';  % assign legend string
h15.DisplayName = 'Blue Tit';  % assign legend string
h16.DisplayName = 'Pigeon';  % assign legend string
h17.DisplayName = 'Greenfinch';  % assign legend string
h18.DisplayName = 'Blackbird';  % assign legend string
h19.DisplayName = 'Robin';  % assign legend string
h20.DisplayName = 'Sparrow';  % assign legend string
h21.DisplayName = 'Blue Tit';  % assign legend string
h22.DisplayName = 'Pigeon';  % assign legend string
h23.DisplayName = 'Greenfinch';  % assign legend string
h24.DisplayName = 'Blackbird';  % assign legend string

h29.DisplayName = 'AMBISONICS';  % assign legend string
h30.DisplayName = 'LOUDSPEAKER (LOW PASSED)';  % assign legend string
h31.DisplayName = 'LOUDSPEAKER (NO REVERB)';  % assign legend string
h32.DisplayName = 'LOUDSPEAKER (SOLO)';  % assign legend string

% Create and format legend
lg = legend([h29 h30 h31 h32 h1 h7 h13 h19...
    h2 h8 h14 h20 ...
    h3 h9 h15 h21 ...
    h4 h10 h16 h22 ...
    h5 h11 h17 h23 ...
    h6 h12 h18 h24 ...
    h25 h26 h27 h28],'Orientation','Horizontal','NumColumns',6,'FontSize',11);
lg.NumColumns = 4;
lg.Location = 'southoutside';

hold off


%-------------------------------------------------------------------------
%% 7. Mann-Whitney U test comparing confidence between ambi and lspk added birds

[pAddedBirdConf,hAddedBirdConf,statsAddedBirdConf] = ranksum(addedBirdsBNambi.Confidence, addedBirdsBNlspk.Confidence)


%-------------------------------------------------------------------------
%% 7. Mann-Whitney U test comparing HB error between ambi and lspk added birds

[pAddedBirdErr,hAddedBirdErr,statsAddedBirdErr] = ranksum(ambiLocAdjErrors(:,3), lspkLPLocAdjErrors(:,3))


%-------------------------------------------------------------------------