% Simulate temporal production EEG
% Apply GLM and regular ERP analysis
%
% Other m-files required:squared_exponential,
% try_gaussian_processed, lsmr, make_stretched_basis,
% uf_vif

% Author: Cameron Hassall, Department of Psychiatry, University of Oxford
% email address: cameron.hassall@psych.ox.ac.uk
% Website: http://www.cameronhassall.com
% March 2020; Last revision: 11-Aug-2022

close all; clear all;

% Set results folder (change as needed)
resultsFolder = '/Users/chassall/OneDrive - Nexus365/Projects/2020_EEG_TimeEstimationTS_Hassall/results';

% Choose simulation
% 'standard' main simulation of manuscript, three durations, two fixed and
% one scaled component
% 'cnv': rather than a scaled component, use a cnv-like response with
% constant slope

whichSimulation = 'standard';

switch whichSimulation
    case {'standard','cnv','flat','pulse','narrowpulse','boxcar','trf','trfparam','trfellapsedt'}
        intervals = [1000 2000 3000];
    case 'production'
        intervals = [800 1650 2500];
    case 'prediction'
        intervals = [700 1300];
end

%% simulate data with two temporally scaled evoked potentials

%non-temporally scaled potentials - these are classic stimulus-locked evoked
%potentials, driven by each cue. It lasts 400ms
stim_evoked = try_gaussian_processes([0:50:400]',[0 1.2 -1.2 0.6 -0.4 0.7 0.5 0.1 0]',[0:400]',25); %stimulus locked potential, evoked to each cue
button_press= try_gaussian_processes([-400:100:100]',[0 0.3 0.8 1.7 2.5 0]',[-400:100]',50); %response locked potential, evoked to each cue

cnv_xs = 0:3001;
cnv_ys = linspace(0,-10,length(cnv_xs));
trf = try_gaussian_processes([0:20:80]',0.03*[0 1.5 3 1.5 0]',[0:80]',25);
switch whichSimulation
    case 'standard'
        cnv_ys_temp = {};
        for i = 1:3
            cnv_ys_temp{i} = try_gaussian_processes([0:intervals(i)/5:intervals(i)]',[0 0.3 0.5 0.3 -0.2 0]',[0:intervals(i)]',intervals(i)/10);
        end
    case 'cnv'
        cnv_ys_temp = {};
        maxCNV = 6;
        cnvSlope = maxCNV/3500;
        cnv_ys_temp{1} = cnvSlope*(0:1500)';
        cnv_ys_temp{2} = cnvSlope*(0:2500)';
        cnv_ys_temp{3} = cnvSlope*(0:3000)';
    case 'flat'
        cnv_ys_temp = {};
        cnv_ys_temp{1} = 3*ones(1501,1);
        cnv_ys_temp{2} = 3*ones(2501,1);
        cnv_ys_temp{3} = 3*ones(3501,1);
    case 'pulse'
        cnv_ys_temp = {};
        cnv_ys_temp{1} = try_gaussian_processes([0:250:1000]',3*[0 1 1 1 0]',[0:1000]',500);
        cnv_ys_temp{2} = try_gaussian_processes([0:250:2000]',3*[0 1 1 1 1 1 1 1 0]',[0:2000]',500);
        cnv_ys_temp{3} = try_gaussian_processes([0:250:3000]',3*[0 1 1 1 1 1 1 1 1 1 1 1 0]',[0:3000]',500);
        plot(cnv_ys_temp{1}); hold on;
        plot(cnv_ys_temp{2});
        plot(cnv_ys_temp{3});
    case 'narrowpulse'
        cnv_ys_temp = {};
        cnv_ys_temp{1} = try_gaussian_processes([0:250:1000]',3*[0 1 1 1 0]',[0:1000]',100);
        cnv_ys_temp{2} = try_gaussian_processes([0:250:2000]',3*[0 1 1 1 1 1 1 1 0]',[0:2000]',100);
        cnv_ys_temp{3} = try_gaussian_processes([0:250:3000]',3*[0 1 1 1 1 1 1 1 1 1 1 1 0]',[0:3000]',100);
        plot(cnv_ys_temp{1}); hold on;
        plot(cnv_ys_temp{2});
        plot(cnv_ys_temp{3});
    case 'boxcar'
        cnv_ys_temp = {};
        cnv_ys_temp{1} = [zeros(1,200) 3*ones(1,601) zeros(1,200)]';
        cnv_ys_temp{2} = [zeros(1,200) 3*ones(1,1601) zeros(1,200)]';
        cnv_ys_temp{3} = [zeros(1,200) 3*ones(1,2601) zeros(1,200)]';
        plot(cnv_ys_temp{1}); hold on;
        plot(cnv_ys_temp{2});
        plot(cnv_ys_temp{3});
    case {'trf','trfparam','trfellapsedt'}
        boxcars = {};
        boxcars{1} = [zeros(1,200) ones(1,601) zeros(1,200)]';
        boxcars{2} = [zeros(1,200) ones(1,1601) zeros(1,200)]';
        boxcars{3} = [zeros(1,200) ones(1,2601) zeros(1,200)]';

        % This might not be efficient, but hopefully it's clear
        cnv_ys_temp = {};
        trfXs = {};
        for s = 1:3
            trfXs{s} = zeros(length(boxcars{s}),length(trf));
            for i = 1:length(boxcars{s})
                if boxcars{s}(i)
                    for j = 1:length(trf)
                        if strcmp(whichSimulation,'trfparam')
                            trfXs{s}((i+(j-1)),j) = s;
                        elseif strcmp(whichSimulation,'trfellapsedt')
                            trfXs{s}((i+(j-1)),j) = i/3001;
                        else
                            trfXs{s}((i+(j-1)),j) = 1;
                        end
                    end
                end
            end
            cnv_ys_temp{s} = trfXs{s}*trf;
        end
        % cspy(thisX);

        
end

% plot(cnv_ys_temp{1} ); hold on;
% plot(cnv_ys_temp{2} ); 
% plot(cnv_ys_temp{3} ); 

%scaled potentials - this is the exact same evoked potential, but it lasts 1000ms
%in response to cue 1 and 2000ms in response to cue 2
tscs = {};
for i = 1:length(intervals)
    tscs{i} = cnv_ys_temp{i}(1:(intervals(i)+1));
end

%% set trial start times
%tstart1 and tstart2 will be lists of when each trial started, in ms

nTrials = 50; %number of trials per condition

for i = 1:length(intervals)
    ICIs(i) = intervals(i) + 2000;
end
blockTime = sum(ICIs);

tstarts = [];

tstarts(1,1) = ICIs(1);
for i = 2:length(intervals)
    tstarts(1,i) = tstarts(1,i-1)+ICIs(i-1);
end

for t = 2:nTrials
    for i = 1:length(intervals)
        tstarts(t,i) = tstarts(t-1,i) + blockTime;
    end
end

% how long should we 'record' simulated data for
recording_length = max(max(tstarts))+max(ICIs);

%% simulated EEG data

%simulate some noise
EEGnoise = 0.5*randn(recording_length,1);

%now add in some signal:
EEGdat = EEGnoise;
EEGdatScaled = EEGdat;
for t = 1:nTrials
    for i = 1:length(intervals)
        EEGdatScaled(tstarts(t,i):tstarts(t,i)+intervals(i)) = EEGdatScaled(tstarts(t,i):tstarts(t,i)+intervals(i)) + tscs{i};
        EEGdat(tstarts(t,i):tstarts(t,i)+400) = EEGdat(tstarts(t,i):tstarts(t,i)+400) + stim_evoked;
        EEGdat(tstarts(t,i):tstarts(t,i)+intervals(i)) = EEGdat(tstarts(t,i):tstarts(t,i)+intervals(i)) + tscs{i};
        EEGdat(tstarts(t,i)+intervals(i)-400:tstarts(t,i)+intervals(i)+100) = EEGdat(tstarts(t,i)+intervals(i)-400:tstarts(t,i)+intervals(i)+100) + button_press;
    end
end

%% first try epoching the data and condition-averaging - "classical" ERP approach

epoched_dat = [];
epochedDatScaled = [];
for t = 1:nTrials
    for i = 1:length(intervals)
        epoched_dat(i,t,:) = EEGdat(tstarts(t,i)-500:tstarts(t,i)+3500);
        epochedDatScaled(i,t,:) = EEGdatScaled(tstarts(t,i)-500:tstarts(t,i)+3500);
    end
end

%% now try convolution modelling instead

r1 = sparse(recording_length,1000); % Stim-locked
r2 = sparse(recording_length,1000); % Response-locked
r3 = sparse(recording_length,1000); % Time-scaled

%regressor r1 captures stim-locked ERP for all cues
for i = 1:1000
    for iInterval = 1:length(intervals)
        r1(tstarts(:,iInterval)-1+i,i) = 1;
    end
end 

%regressor r2 captures resp-locked ERP for all cues
for i = 1:1000
    for iInterval = 1:length(intervals)
        r2(tstarts(:,iInterval)+intervals(iInterval)-800+i,i) = 1;
    end
end 

%regressor r3 captures temporally scaled potential for both cues
for j = 1:size(tstarts,1)
    for iInterval = 1:length(intervals)
        r3(tstarts(j,iInterval):tstarts(j,iInterval)+intervals(iInterval)-1,:) = make_stretched_basis(intervals(iInterval),1000);
    end
end 

%also consider r4, a fixed-time signal that starts at 0 and changes the
%same way regardless of the duration
r4 = sparse(recording_length,3000);
for j = 1:size(tstarts,1)
    for iInterval = 1:length(intervals)
        r4(tstarts(j,iInterval):tstarts(j,iInterval)+intervals(iInterval)-1,1:intervals(iInterval)) = eye(intervals(iInterval));
    end
end 

% also consider the method of TRFs, assuming we know when the TRF is active
if strcmp(whichSimulation,'trf') || strcmp(whichSimulation,'trfparam')
    r5 = sparse(recording_length,length(trf));
    for j = 1:size(tstarts,1)
        for iInterval = 1:length(intervals)
            r5(tstarts(j,iInterval):tstarts(j,iInterval)+intervals(iInterval),1:length(trf)) = trfXs{iInterval};
        end
    end
    Xtrf = [r1 r2 r5];
    lsmriterations = 400;
    [betaTRF,ISTOP,ITN] = lsmr(Xtrf,EEGdat,[],10^-8,10^-8,[],lsmriterations);
end


%design matrix
X = [r1 r2 r3];

thisVIF =uf_vif(X);

%betas from regression model:
% See UNFOLD toolbox
lsmriterations = 400;
[beta,ISTOP,ITN] = lsmr(X,EEGdat,[],10^-8,10^-8,[],lsmriterations);

% With only fixed components
XFixedOnly = [r1 r2 r4];
[betaFixedOnly,ISTOP,ITN] = lsmr(XFixedOnly,EEGdat,[],10^-8,10^-8,[],lsmriterations);

% pDM = pinv_reg(X,100000,'onediff',[1000 2000]);
% beta = pDM * EEGdat;

XNoTS = [r1 r2];
[betaNoTS,ISTOP,ITN] = lsmr(XNoTS,EEGdat,[],10^-8,10^-8,[],lsmriterations);

residuals = EEGdat - X * beta;
residualsFixedOnly = EEGdat - XFixedOnly * betaFixedOnly;

errorTS = mean(residuals.^2);
errorFixedOnly = mean(residualsFixedOnly.^2);

% residualsNoTS = EEGdat - XNoTS * betaNoTS;
% residualsTest = EEGdat - XNoTS * beta(1:2000);

save(fullfile(resultsFolder,['task1-2_simulationresults_' whichSimulation '.mat']));