%% marshmallow_header
% Robert Molitor, January 2020
% header file for the experiment testing whether "expertise" improves mnemonic disctrimination
% creates a subject data folder, mat file, and workspace structure called 'header' that is read into sub functions
%
% subfunctions: marshmallow_run, marshmallow_studytest
%
% ##input##
% none immediately, but will prompt:
% subject number
% subject age
% subject sex
%
% ##output##
% data directory, data/marshmallow_xx_yy
% header file, data/marshmallow_xx_yy/marshmallow_xx_yy_header.mat
% workspace structure, header, that contains all of the experiment parameters for subfunctions
%
% ##current design (as of 01/27/20)##
% 2 scene categories (beaches and gazebos)
% 80 scenes per category
% nback with one of the categories, associative learning with both categories
% post-test including scene-object pair learning and interference test
%
% --nback run--
% 2back task
% 20 runs
% 50 trials per run, 10 items are 'yes' nback responses
% 2 repetitions per item per run
% 1s stimulus duration
% 1s ITI
% run length: ~110s (~36m total)
%
% --pairmate selection for post-test--
% 20 pairs per category (40 total), randomly selected
%
% --study and test--
% 3 study and 3 test runs
% 80 trials per phase
% study: 1.5s stimulus x 2 + 0.5s ISI + 2s ITI (7m 20s per study run)
% test: 1.5s stimulus + 0.5s ISI + 5s memory array (target, competitor, non-competitor) + 2s ITI (12m max per test run)

%% clear any preexisting stuff in the workspace and command window
clear all %#ok<*CLALL>
clc

%add path to psychtoolbox
addpath(genpath('/Users/kuhllab/Desktop/Psychtoolbox'));

%% SETUP
expname = 'marshmallow';
expversion = '1.0';
exptype = 'behavior_pilot';
header = struct('exp', expname, 'version', expversion, 'exptype', exptype, 'parameters', struct); %experiment info
header.path.exp = pwd; %path for experiment files
par = struct;
header.path.data = [header.path.exp '/data/']; %path for data folder
header.path.stim = [header.path.exp '/stimuli/']; %path for stimuli
header.setuptime = fix(clock); %time this was run

%% SUBJECT INFO
header.subjnum = input('Participant number:  ','s'); %subject number

%make sure data for this subject # doesn't already exist
header.subjinfo = [expname '_',header.subjnum]; %add subject info to header
header.path.subjinfo = [header.path.data sprintf('%s', header.subjinfo)];

if exist(header.path.subjinfo,'file') %abort experiment, rather than overwrite existing data
    disp('Data for this participant already exists. Ending...');
    clear all
    return
end

%figure out the counterbalancing (based on whether the subject # is odd or even)
if rem(str2double(header.subjnum),2) ~= 0
    header.cb = 1;
    header.design.nback.category = 1;
elseif rem(str2double(header.subjnum),2) == 0
    header.cb = 2;
    header.design.nback.category = 2;
end

%if everything is good, make the data folder
mkdir(sprintf('%s/%s',header.path.data, header.subjinfo));

%seed the random number generator by subject number
rng(str2double(header.subjnum));

%% THINGS YOU MIGHT WANT TO EDIT
% PSYCHTOOLBOX STUFF
par.ptb.txtcol = [255 255 255]; %black text
par.ptb.bgcol = [128 128 128]; %white background
par.ptb.txtsize_fix = 75; %size of the fixation cross
par.ptb.txtsize_inst = 50; %size of the text for instructions
par.ptb.txtsize = 75; %text size for almost everything
par.ptb.buffer_test = 50; %buffer between items during the 3AFC test
par.ptb.txtsize_test = 30; %size of test instructions

% TRIALS AND RUNS AND STIM
% RUNS, TRIALS, REPS
%general
par.design.ncategories = 2; %number of scene categories
par.design.nstim = 80; %number of scenes per category

%study and test
par.design.studytest.ncond = 2; %number of conditions
par.design.studytest.nruns = 3; %number of runs
par.design.studytest.nruns_prac = 1; %number of runs in the practice
par.design.studytest.npairs_max = par.design.nstim/2; %maximum number of pairs per category
par.design.studytest.npairs_prcnt = 0.25; %percentage of max pairs to use (provides some cutoff to better distinguish conditions if <1) %PK edit: change to .25
par.design.studytest.npairs_prcnt_bycat = 5; %pairs per condition per category using the cutoff %PK edit: this math wont work with.25, so I have to hardcode
par.design.studytest.ntrls_study = 40; %number of study trials
par.design.studytest.ntrls_test = par.design.studytest.ntrls_study; %number of test trials

%for study and test practice, just hard code the trials
%uses three practice-only scenes (unrelated to the categories in the main task) and three practice-only objects
par.studytest.study_prac = [1 1; 2 2; 3 3]; %first column is practice scene, second column is practice object
par.studytest.test_prac = [1 2 1 3; 2 3 2 1; 3 1 2 3]; %first column is cue, columns two through four are the objects

%nback
par.design.nback.nback = 2; %how many back for the nback
par.design.nback.nback_prcnt = 0.2; %percentage of trials that are yes responses
par.design.nback.item_trls = par.design.nstim/2; %number of non nback trials per run, equal to half of all the stimuli in a category
par.design.nback.total_trls = par.design.nback.item_trls/(1 - par.design.nback.nback_prcnt); %total # of trials in an nback run, including the nback
par.design.nback.nback_trls = par.design.nback.total_trls - par.design.nback.item_trls; %actual # of nback trials
par.design.nback.nruns = par.design.nback.item_trls/par.design.nback.nback_trls; %how many runs are needed to do all of the stimuli as the repeats
par.design.nback.ncycles = 5; %how many times to repeat the number of runs necessary to do all of the stimuli as repeats

% TIMING (s)
%exposure
par.timing.nback.stim = 1; %duration of nback items
par.timing.nback.ITI = 1; %iti

%study and test
par.timing.studytest.study_stim = 1.5; %duration of items in the study phase
par.timing.studytest.study_ISI = 0.5; %time between items
par.timing.studytest.study_ITI = 2; %time between study trials
par.timing.studytest.test_cue = 1.5; %duration of cue item
par.timing.studytest.test_ISI = 0.5; %time between cue and selection array
par.timing.studytest.test_array = 5; %time (max) of test array
par.timing.studytest.test_ITI = 2;  %time between test trials

%nback
par.timing.nback.stim = 1; %duration of stimulus items + response deadline
par.timing.nback.ITI = 1; %time between trials

% STIMULI
%image sizing
par.stim.objheight = 256; %object height for movies in pixels
par.stim.objwidth = 256; %object width for movies in pixels

%labels
par.stim.labels = {'beach' 'gazebo'};

%objects
par.stim.objects = randsample(80,40,'false');

%% PAIR LIST
%make pairmates (decided randomly)

%for making random pairs, create all of the possible combinations
combos = nchoosek(1:par.design.nstim,2);
combos = [transpose(1:size(combos,1)) combos]; %give each pair and ID for easier stuff later

%now choose the appropriate number of pairs for each category
pairmates = cell(1,par.design.ncategories);
nbackpairmates = cell(1,par.design.ncategories);
for c = 1:par.design.ncategories
    
    %start the process by randomizing the pairs
    randomizer = Shuffle(1:size(combos,1));
    combos_cat = combos(randomizer,:);
    
    %then get all the pairs needed
    for p = 1:20 %here?
        
        %take the top item as the next pair
        thispair = combos_cat(1,:);
        pairmates{c} = [pairmates{c}; thispair];
        
        %remove other possible pairs that contain those items
        check1 = combos_cat(:,2:3) == thispair(2); %pairs that have the first item
        check2 = combos_cat(:,2:3) == thispair(3); %pairs that have the second item
        checks = sum([check1 check2],2); %pairs that contain any of the itemschecks
        combos_cat = combos_cat(checks == 0,:); %remove pairs
        
    end %end p
    
    %add in condition number, which for random pairs = 0
    %also add in category number
    pairmates{c} = [pairmates{c} zeros(size(pairmates{c},1),1) repmat(c,size(pairmates{c},1),1)];
    nbackpairmates{c} = pairmates{c};
    tmp = pairmates{c};
    pairmates{c} = tmp(randsample(20,10),:);
end %end c

%save it
par.studytest.pairmates = pairmates;
par.nback.pairmates = nbackpairmates;
%% NBACK LIST
%based on the pairs used in the memory test, get the items for the nback
studytest_items = nbackpairmates{header.design.nback.category}; %pairs for the category used in the nback
studytest_items = [studytest_items(:,2); studytest_items(:,3)]; %get all items for the category

%remove items from the memory test from the nback
nback_items = transpose(1:par.design.nstim);
for nb = 1:size(studytest_items,1)
    
    nback_items(nback_items == studytest_items(nb)) = [];
    
end %end nb

%% ITEM LIST
%abstracted list of items that contains individual pair items and objects
itemlist = []; %for making the study list
itemlist_bycat = cell(1,par.design.ncategories); %for the test (slightly easier indexing with this setup)
objects = reshape(Shuffle(par.stim.objects),[],par.design.ncategories); %randomly assign objects to pair items
for c = 1:par.design.ncategories
    
    for p = 1:par.design.studytest.npairs_prcnt_bycat * par.design.studytest.ncond %here?
        
        itemlist = [itemlist; repmat([c p],2,1) transpose([1 2])]; %#ok<*AGROW> %category, pair #, item #
        
    end %end p
    
    %separate out by category and add in objects
    itemlist_bycat{c} = [itemlist(itemlist(:,1) == c,:) objects(:,c)];
    
end %end c

%add in objects to the whole list
itemlist = [itemlist objects(:)];

%add in competitor objects for ease later
cats = 1:par.design.ncategories;
competitors = zeros(size(itemlist,1),1);
for t = 1:size(itemlist,1)
    
    %really complicated way to just find the other pairmate
    %basically, find item in the same category that is the other pairmate to this pair
    index = find(itemlist(:,1) == itemlist(t,1) & itemlist(:,3) == cats(cats ~= itemlist(t,3)) & itemlist(:,2) == itemlist(t,2));
    
    %get the competitor for this pair and add it to the item list
    competitors(t) = itemlist(index,4);
    
end

%add competitors to the item list
itemlist = [itemlist competitors];
competitors_bycat = reshape(competitors,[],par.design.ncategories);
for c = 1:par.design.ncategories
    
    %separate out by category and add in objects
    itemlist_bycat{c} = [itemlist_bycat{c} competitors_bycat(:,c)];
    
end %end c

%% STUDY
%create trial onsets
%randomize order of study trials
%maximum number of same category trials in a row is 4
%pseudoradomize so there are at least 2 intervening items between pairmates

%onsets first
onsets = cell(1,par.design.studytest.nruns);
for r = 1:par.design.studytest.nruns
    
    thisonset = 2; %start after a short lag
    
    for t = 1:(par.design.studytest.ntrls_study-1) %here?
        
        thisonset = [thisonset; thisonset(end) + 2*par.timing.studytest.study_stim + par.timing.studytest.study_ISI + par.timing.studytest.study_ITI]; %time between study trials]; %based on previous onset + two stimuli presentations + ISI + jitter
        
    end %end t
    
    onsets{r} = thisonset;
    
end %end r

%save out some onsets for the practice
par.studytest.onsets_prac = onsets{1}(1:size(par.studytest.study_prac,1));

%now trial order
runs = cell(1,par.design.studytest.nruns);

for r = 1:par.design.studytest.nruns
    
    repcheck = 1;
    while repcheck
        
        repcheck = 0; %reset repcheck every iteration of the while loop
        
        %randomize trial order
        randomizer = Shuffle(1:size(itemlist,1));
        temp = itemlist(randomizer,:);
        
        %check for repeats of category items
        for t = 5:size(temp,1)
            
            if temp(t,1) == temp(t-1,1) && temp(t,1) == temp(t-2,1) && temp(t,1) == temp(t-3,1) && temp(t,1) == temp(t-4,1)
                repcheck = repcheck + 1;
            end
            
        end %end t
        
        %now check for proximity of pairmates of the same category
        for t = 3:size(temp,1)
            
            if temp(t,2) == temp(t-1,2) || temp(t,2) == temp(t-2,2) %pairmate check
                if temp(t,1) == temp(t-1,1) || temp(t,1) == temp(t-2,1) %category check
                    repcheck = repcheck + 1;
                end
            end
            
        end %end t
        
    end %end while
    
    %add on onsets
    runs{r} = [temp onsets{r}];
    
end %end r

%save it
par.studytest.runs_study = runs; %final study matrix: category, pair #, item #, target, competitor, onset

%% TEST
%randomize order of test trials
%pseudoradomize so there are at least 2 intervening items between pairmates
%maximum number of same category trials in a row is 4
%figure out the third non-competitor item taken from the other category, which does not change across runs
%counterbalance positioning

%start with the positioning
pos = 1:3; %1 = target, 2 = competitor, 3 = non-competitor; positions left, middle, right
all_pos = perms(pos); %all permutations of item positioning
pos_rep = par.design.studytest.ntrls_test*par.design.studytest.nruns/size(all_pos,1); %how many times the permutations can be repeated across the test

%check to see if they # of trials are divisible by the number of positions for an equal distribution (checks for whole number through rounding)
if ceil(pos_rep) ~= pos_rep %if not, round up to have extra positions available
    pos_test = repmat(all_pos,ceil(pos_rep),1); %permuations sampled across the entire test
else
    pos_test = repmat(all_pos,pos_rep,1); %permuations sampled across the entire test
end

%evenly distribute across runs
repcheck = 1;
while repcheck
    
    repcheck = 0; %restart each iteration at zero
    
    randomizer = Shuffle(1:size(pos_test,1));
    pos_rand = pos_test(randomizer,:); %shuffled positions
    
    %split positions by run
    pos_byrun = [];
    for p = 1:length(pos)
        pos_byrun = cat(2,pos_byrun,reshape(pos_rand(:,p),par.design.studytest.ntrls_test,1,par.design.studytest.nruns));
    end
    
    % pos_byrun = reshape(pos_rand,par.design.studytest.ntrls_test,length(pos),par.design.studytest.nruns); %positions split up by run
    
    poscheck = zeros(par.design.studytest.nruns,length(pos));
    for r = 1:par.design.studytest.nruns
        for p = 1:length(pos)
            poscheck(r,p) = sum(pos_byrun(:,p,r) == 1); %how many times the correct answer is in each position each run
        end
    end
    
    %try to keep the imbalance of correct answer positioning minimal across runs
    if max(max(poscheck,[],2)) > (floor(par.design.studytest.ntrls_test/length(pos)) + 1) || min(min(poscheck,[],2)) < (floor(par.design.studytest.ntrls_test/length(pos)) - 1)
        repcheck = repcheck + 1;
    end
    
    
end %end repcheck

%now assign non-competitor items to each item
%should work if there is more than two categories
repcheck = 1;
while repcheck
    
    repcheck = 0;
    
    %randomly shuffle objects to new categories
    cats = 1:par.design.ncategories;
    randomizer = Shuffle(cats);
    
    %if any object list has the same category label, re-shuffle
    if sum(cats == randomizer) > 0
        repcheck = repcheck + 1;
    end
    
end %end repcheck

%assign objects to new categories and randomize
%shuffle does each column independently, so only one call is needed here
objshuff = objects(:,randomizer);
objshuff = Shuffle(objshuff);

%add in non-competitors to the item list
%itemlist is now: category, pair #, item #, target, competitor, non-competitor object
itemlist = [itemlist objshuff(:)];
for c = 1:par.design.ncategories
    
    itemlist_bycat{c} = [itemlist_bycat{c} objshuff(:,c)];
    
end %end c

%randomize trial orders
runs = cell(1,par.design.studytest.nruns);

for r = 1:par.design.studytest.nruns
    
    repcheck = 1;
    while repcheck
        
        repcheck = 0; %reset repcheck every iteration of the while loop
        
        %randomize trial order
        randomizer = Shuffle(1:size(itemlist,1));
        temp = itemlist(randomizer,:);
        
        %check for repeats of category items
        for t = 5:size(temp,1)
            
            if temp(t,1) == temp(t-1,1) && temp(t,1) == temp(t-2,1) && temp(t,1) == temp(t-3,1) && temp(t,1) == temp(t-4,1)
                repcheck = repcheck + 1;
            end
            
        end %end t
        
        %now check for proximity of pairmates of the same category
        for t = 3:size(temp,1)
            
            if temp(t,2) == temp(t-1,2) || temp(t,2) == temp(t-2,2) %pairmate check
                if temp(t,1) == temp(t-1,1) || temp(t,1) == temp(t-2,1) %category check
                    repcheck = repcheck + 1;
                end
            end
            
        end %end t
        
        %finally, make sure there are not trials with the same objects back to back
        for t = 3:size(temp,1)
            
            if sum(ismember(temp(t,4:6),temp(t-1,4:6))) > 0
                repcheck = repcheck + 1;
            end
            
        end %end t
        
    end %end while
    
    %add in position information and save it
    runs{r} = [temp pos_byrun(:,:,r)]; %final test matrix: category, pair #, item #, target, competitor, non-competitor, left item, middle item, right item
    
end %end r

%save it
par.studytest.runs_test = runs;

%% N-BACK
%make all of the runs; first column is the stimulus, second column is the correct response (1 = no repeat, 2 = repeat)

%initialize empty matrix to hold all of the run data
runs = [];

%for each cycle, make enough runs so that every item can serve as an nback repeat
for c = 1:par.design.nback.ncycles
    
    nback_byrun = reshape(Shuffle(nback_items),[],par.design.nback.nruns); %figure out which stimuli will be the nback items in each run
    
    for r = 1:par.design.nback.nruns
        
        %make the nback list, check to make sure the same item isn't repeated in the nback
        repcheck = 1;
        
        while repcheck
            
            repcheck = 0; %reset repcheck every iteration
            
            run_start = ones(par.design.nback.nback,1); %first few trials that cannot actually be nback
            nback_list = vertcat(ones(par.design.nback.item_trls - par.design.nback.nback,1),ones(par.design.nback.nback_trls,1)*2); %add in the trials beyond the requisite starters as well as the nback trials for randomizing
            trl_code = vertcat(run_start,Shuffle(nback_list)); %bring it together
            
            for t = par.design.nback.nback + 1:size(trl_code,1)
                
                if trl_code(t) == trl_code(t - par.design.nback.nback) && trl_code(t) == 2 %checks for nback yes responses for the same items repeatedly
                    repcheck = repcheck + 1;
                end
                
            end %end t
            
        end %end while repcheck
        
        %get the nback items for this run and create what will be the final list
        thisnback = nback_byrun(:,r);
        trlstim = nan(par.design.nback.total_trls,1);
        
        %remove the nback items from the rest of the list
        thislist = Shuffle(nback_items);
        
        for nb = 1:size(thisnback,1)
            
            thislist(thislist == thisnback(nb)) = [];
            
        end %end nb
        
        %add in the nback trials first
        for t = 1:par.design.nback.total_trls
            
            if trl_code(t) == 2
                
                trlstim(t) = thisnback(1);
                trlstim(t - par.design.nback.nback) = thisnback(1);
                thisnback(1) = [];
                
            end
            
        end %end t
        
        %add in the other items
        for t = 1:par.design.nback.total_trls
            
            if isnan(trlstim(t))
                trlstim(t) = thislist(1);
                thislist(1) = [];
            end
            
        end %end t
        
        %bring it together
        thisrun = [trlstim trl_code];
        
        runs = cat(3,runs,thisrun);
        
    end %end r runs
    
end %end c cycles

%save it
par.nback.runs = runs;

%make the practice
par.nback.runs_prac = [transpose([4 5 4 6 7 6 8 9 10 11 12 11 13]) transpose([1 1 2 1 1 2 1 1 1 1 1 2 1])];

%finally, make the onsets, which will be the same across all runs
thisonset = 2; %start after a short lag

for t = 1:(par.design.nback.total_trls - 1)
    
    thisonset = [thisonset; thisonset(end) + par.timing.nback.stim + par.timing.nback.ITI];
    
end %end t

%save it
par.nback.onsets = thisonset;

%% SAVE IT ALL
header.parameters = par;
clearvars -except header;
save(sprintf('%s/%s_header', header.path.subjinfo, header.subjinfo))

%fin