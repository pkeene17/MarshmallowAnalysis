function data = marshmallow_studytest(header,exp)
%% marshmallow_studytest
% Robert Molitor, January 2020
% Edited by Paul Keene, October 2021
% study/test for marshmallow2-electric_marshmalloogaloo
%
% ##input##
% header: header created from marshmallow_header
% exp: 0 = practice, 1 = main experiment
% EXAMPLE: marshmallow_studytest(header,1) will run the main task
%
% ##study output##
% run
% trial
% cb: 1 = beaches nback, 2 = gazebos nback
% category: 1 = beach, 2 = gazebo
% pairnum: which pair within a category this item comes from
% pairitem: which item is being shown on this trial (pairmate 1 or pairmate 2)
% pairid: ID # of this pair (based on all possible combinations)
% scenenum: which scene is being shown
% objnum: which object is being shown
% condition: condition (0-4)
% onset: trial onset
%
% ##test output##
% run
% trial
% cb: 1 = beaches nback, 2 = gazebos nback
% category_cue: category of the cued scene and competitor (1 = beach, 2 = gazebo)
% category_noncomp: category of the non-competitor item (1 = beach, 2 = gazebo)
% pairnum_cue: which pair within a category the cue and competitor comes from
% pairnum_noncomp: which pair within a category the non-competitor comes from
% pairitem_cue: which item is being shown on this trial (pairmate 1 or pairmate 2) for the cue
% pairitem_noncomp: which item is being shown on this trial (pairmate 1 or pairmate 2) for the non-competitor
% pairid_cue: ID # of this pair (based on all possible combinations) for the cue and competitor
% pairid_noncomp: ID # of this pair (based on all possible combinations) for the non-competitor
% scenenum: which scene is being shown for the cue
% objnum_obj1: which object is being shown for the left item
% objnum_obj2 which object is being shown for the middle item
% objnum_obj3: which object is being shown for the right item
% condition_cue: similarity condition (0-4) of the cue and competitor
% condition_noncomp: similarity condition (0-4) of the non-competitor
% objpos1: which item appeared on the left (1 = target, 2 = competitor, 3 = non-competitor)
% objpos2: which item appeared on the middle (1 = target, 2 = competitor, 3 = non-competitor)
% objpos3:which item appeared on the right (1 = target, 2 = competitor, 3 = non-competitor)
% resp: actual response
% score: 1 = target, 2 = competitor, 3 = non-competitor
% rt: response time (s)
%

%% READY, SET, GO
par = header.parameters;
data = struct('subjnum', header.subjnum);
data.header = header;
addpath(header.path.stim);
addpath(genpath('/Users/kuhllab/Desktop/Psychtoolbox'));

%% OUTPUT TEXT FILES
if exp %only make textfile outputs for the main task
    %study file
    outSname = sprintf('%s/%s_study_%d',header.path.subjinfo,header.subjinfo,exp);
    if exist([outSname '.txt'],'file') == 2
        error('The data file for this task already exists.  Please check your input parameters.');
    end
    
    Sid=fopen([outSname '.txt'], 'w'); % open the file
    fprintf(Sid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
        'run','trial','cb','category','pairnum','pairitem','pairid','scenenum','objnum','condition','onset'); %create header
    
    %test file
    outTname = sprintf('%s/%s_test_%d',header.path.subjinfo,header.subjinfo,exp);
    if exist([outTname '.txt'],'file') == 2
        error('The data file for this task already exists.  Please check your input parameters.');
    end
    
    Tid=fopen([outTname '.txt'], 'w'); % open the file
    fprintf(Tid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
        'run','trial','cb','category_cue','category_noncomp','pairnum_cue','pairnum_noncomp','pairitem_cue','pairitem_noncomp','pairid_cue','pairid_noncomp','scenenum','objnum_obj1','objnum_obj2','objnum_obj3','condition_cue','condition_noncomp','objpos1','objpos2','objpos3','resp','score','rt'); %create header
end

%% SCREEN SETUP
oldEnableFlag = Screen('Preference','SuppressAllWarnings');
Screen('Preference', 'VisualDebugLevel', 1);
Screen('Preference', 'SkipSyncTests', 2); %will not show warning screen
allscreens = Screen('screens'); %how many screens are connected to the computer
screennum = max(allscreens); %choose the last one (not the main screen)
[window, screenrect] = Screen('OpenWindow',screennum,par.ptb.bgcol);
Screen('BlendFunction',window,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA); %allows for the transparent backgrounds on the PNG files
screenLeft = screenrect(RectLeft); %pixel value for left of screen
screenRight = screenrect(RectRight); %pixel value for right of screen
screenTop = screenrect(RectTop); %pixel value for top of screen
screenBottom = screenrect(RectBottom); %pixel value for bottom of screen
screenWidth = screenRight - screenLeft; %pixel value for screen width
ycenter = (screenBottom - screenTop)/2; %ycenter pixel value
xcenter = (screenRight - screenLeft)/2; %xcenter pixel value
Screen(window, 'TextSize', par.ptb.txtsize);
Screen(window, 'TextFont', 'Arial');
HideCursor;

%% IMAGE SIZING
objheight = par.stim.objheight; %height in pixels
objwidth = par.stim.objwidth; %width in pixels

%% STUDY POSITIONING (ALSO TEST CUE POSITIONING)
%image displayed in the center of the screen
Sx1 = xcenter - objwidth/2; %left side of image
Sy1 = ycenter - objheight/2; %top of image
Sx2 = xcenter + objwidth/2; %right side of image
Sy2 = ycenter + objheight/2; %bottom of image
itemrect = [Sx1 Sy1 Sx2 Sy2]; %rectangular coordinates for the object location

%% TEST POSITIONING
%three items, positioned vertically centered
%object 1
O1x1 = xcenter - objwidth - par.ptb.buffer_test - objwidth/2; %left side of image
O1y1 = ycenter - objheight/2; %top of image
O1x2 = xcenter - par.ptb.buffer_test - objwidth/2; %right side of image
O1y2 = ycenter + objheight/2; %bottom of image
obj1rect = [O1x1 O1y1 O1x2 O1y2]; %rectangular coordinates for the object location

%object 2
O2x1 = screenWidth/2 - objwidth/2; %left side of image
O2y1 = ycenter - objheight/2; %top of image
O2x2 = screenWidth/2 + objwidth/2; %right side of image
O2y2 = ycenter + objheight/2; %bottom of image
obj2rect = [O2x1 O2y1 O2x2 O2y2]; %rectangular coordinates for the object location

%object 3
O3x1 = xcenter + par.ptb.buffer_test + objwidth/2; %left side of image
O3y1 = ycenter - objheight/2; %top of image
O3x2 = xcenter + objwidth + par.ptb.buffer_test + objwidth/2; %right side of image
O3y2 = ycenter + objheight/2; %bottom of image
obj3rect = [O3x1 O3y1 O3x2 O3y2]; %rectangular coordinates for the object location

%% SET RESPONSE OPTIONS
leftkey = KbName('j'); %button to select the left object
middlekey = KbName('k'); %button to select the middle object
rightkey = KbName('l'); %button to select the right object

%% RUN LOOP
if exp
    runs = par.design.studytest.nruns;
else
    runs = par.design.studytest.nruns_prac;
end

for r = 1:runs
    
    %% RUN START
    data.startTime = fix(clock);
    
    %% MAKE EMPTY MATRICES TO HOLD ALL OF THE TRIAL DATA
    studydata = [];
    testdata = [];
    
    %% CHOOSE STIMULI FOR THIS RUN
    run = r;
    if exp %main task
        studylist = par.studytest.runs_study{run}; %study stim
        testlist = par.studytest.runs_test{run}; %test stim
    else %practice
        studylist = par.studytest.study_prac;
        testlist = par.studytest.test_prac;
    end
    
    %% STUDY START SCREEN
    Screen(window, 'FillRect', par.ptb.bgcol, screenrect);
    Screen(window, 'TextSize', par.ptb.txtsize_inst);
    
    if exp %real deal
        studytext = sprintf('STUDY %d/%d\n\nPress J, K, or L to start!',run,runs);
    else %practice
        studytext = sprintf('PRACTICE STUDY %d/%d\n\nPress J, K, or L to start!',run,runs);
    end
    
    DrawFormattedText(window,studytext,'center','center',par.ptb.txtcol);
    Screen(window, 'Flip');
    WaitSecs(1);
    
    %wait for button press to start the experiment
    clear keyCode;
    clear keyIsDown;
    
    respwait = 1;
    while respwait
        [keyIsDown, ~, keyCode] = KbCheck; %check for key response from the participant
        if keyIsDown
            if keyCode(leftkey) || keyCode(middlekey) || keyCode(rightkey)
                break
            end %end of keycode check
        end  %end of kbcheck
    end %end of while
    
    %% STUDY LOOP
    %initial fixation
    Screen(window, 'TextSize', par.ptb.txtsize_fix);
    Screen(window, 'FillRect', par.ptb.bgcol, screenrect);
    DrawFormattedText(window,'+','center','center',par.ptb.txtcol); %central fixation during ITI
    Screen(window, 'Flip');
    
    %start of the study phase
    data.studyStart = fix(clock);
    studystart = GetSecs;
    
    %loop through all the trials
    for t = 1:size(studylist,1)
        
        %figure out which stimuli should be shown and the relevant trial information
        if exp %main experiment
            
            %variables of interest for study trials
            category_cue = studylist(t,1);
            pairnum = studylist(t,2);
            pairitem = studylist(t,3);
            pairid = par.studytest.pairmates{category_cue}(pairnum,1);
            scenenum = par.studytest.pairmates{category_cue}(pairnum, pairitem + 1); %pair item offset by 1 because the first column in pairmates is the pairid
            objnum = studylist(t,4);
            condition_cue = par.studytest.pairmates{category_cue}(pairnum,5);
            onset = studylist(t,6);
            
        else %practice
            
            scenenum = studylist(t,1);
            objnum = studylist(t,2);
            onset = par.studytest.onsets_prac(t);
            
        end %end exp check
        
        %get the label for the scene
        if scenenum < 10
            snum = ['00' num2str(scenenum)];
        elseif scenenum < 100
            snum = ['0' num2str(scenenum)];
        else
            snum = num2str(scenenum);
        end
        
        if exp
            stim1 = [par.stim.labels{category_cue} '_' snum '.jpg'];
        else
            stim1 = ['prac_scene_' snum '.jpg'];
        end
        
        %get the label for the object
        if objnum < 10
            onum = ['00' num2str(objnum)];
        elseif objnum < 100
            onum = ['0' num2str(objnum)];
        else
            onum = num2str(objnum);
        end
        
        if exp
            stim2 = ['object_' onum '.jpg'];
        else
            stim2 = ['prac_object_' onum '.jpg'];
        end
        
        %figure out when to show the first stimulus
        stim1time = studystart + onset;
        
        %load in the stimuli
        [stim1vals, ~, ~] = imread(stim1);
        [stim2vals, ~, ~] = imread(stim2);
        
        %show the first item and figure out when to flip back to fixation
        Screen('PutImage', window, stim1vals, itemrect); %show left object
        flipped = Screen(window, 'Flip',stim1time);
        ISItime = flipped + par.timing.studytest.study_stim;
        
        %flip back to fixation for the ISI and figure out when to show the second stimulus
        Screen(window, 'TextSize', par.ptb.txtsize_fix);
        Screen(window, 'FillRect', par.ptb.bgcol, screenrect);
        DrawFormattedText(window,'+','center','center',par.ptb.txtcol); %central fixation during ITI
        flipped = Screen(window, 'Flip',ISItime);
        stim2time = flipped + par.timing.studytest.study_ISI;
        
        %show the second item figure out when to flip back to fixation
        Screen('PutImage', window, stim2vals, itemrect); %show left object
        flipped = Screen(window, 'Flip',stim2time);
        ITItime = flipped + par.timing.studytest.study_stim;
        
        %flip to fixation for the ITI
        Screen(window, 'TextSize', par.ptb.txtsize_fix);
        Screen(window, 'FillRect', par.ptb.bgcol, screenrect);
        DrawFormattedText(window,'+','center','center',par.ptb.txtcol); %central fixation during ITI
        Screen(window, 'Flip',ITItime);
        
        %save trial into matrix and write to text file
        if exp
            thistrl = [run t header.cb category_cue pairnum pairitem pairid scenenum objnum condition_cue onset];
            studydata = [studydata; thistrl]; %#ok<AGROW>
            fprintf(Sid,'%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.1f\n',thistrl);
        end
        
    end %end of t trials
    
    %wait for the final fixation
    WaitSecs(par.timing.studytest.study_ITI);
    
    %end of study phase
    data.studyStop = fix(clock);
    
    %% BETWEEN RUN INTERVAL
    Screen(window, 'TextSize', par.ptb.txtsize_fix);
    Screen(window, 'FillRect', par.ptb.bgcol, screenrect);
    DrawFormattedText(window,'+','center','center',par.ptb.txtcol);
    Screen(window, 'Flip');
    WaitSecs(2);
    
    %% TEST START SCREEN
    Screen(window, 'TextSize', par.ptb.txtsize_test);
    if exp
        testtext = sprintf('TEST %d/%d\n\nPress J for the left item, K for the middle item, and L for the right item.\n\nPress J, K, or L to start the test!',run,runs);
    else
        testtext = sprintf('PRACTICE TEST %d/%d\n\nPress J for the left item, K for the middle item, and L for the right item.\n\nPress J, K, or L to start the test!',run,runs);
    end
    
    DrawFormattedText(window,testtext,'center','center',par.ptb.txtcol);
    Screen(window, 'Flip');
    Screen(window, 'TextSize', par.ptb.txtsize);
    WaitSecs(1);
    
    %wait for button press to start the experiment
    clear keyCode;
    clear keyIsDown;
    
    respwait = 1;
    while respwait
        [keyIsDown, ~, keyCode] = KbCheck; %check for key response from the participant
        if keyIsDown
            if keyCode(leftkey) || keyCode(middlekey) || keyCode(rightkey)
                break
            end %end of keycode check
        end  %end of kbcheck
    end %end of while
    
    %% TEST LOOP
    %initial fixation
    Screen(window, 'FillRect', par.ptb.bgcol, screenrect);
    DrawFormattedText(window,'+','center','center',par.ptb.txtcol);
    Screen(window, 'Flip');
    WaitSecs(2);
    
    %start of the test phase
    data.testStart = fix(clock);
    
    %go through all the trials
    for t = 1:size(testlist,1)
        
        %figure out which stimuli should be shown and the relevant trial information
        if exp
            
            %for the cue
            category_cue = testlist(t,1);
            pairnum_cue = testlist(t,2);
            pairitem_cue = testlist(t,3);
            pairid_cue = par.studytest.pairmates{category_cue}(pairnum_cue,1);
            scenenum = par.studytest.pairmates{category_cue}(pairnum_cue, pairitem_cue + 1); %pair item offset by 1 because the first column in pairmates is the pairid
            condition_cue = par.studytest.pairmates{category_cue}(pairnum_cue,5);
            objnum = testlist(t,4:6); %target, competitor, non-competitor
            
            %for the non-competitor
            category_noncomp = testlist(testlist(:,4) == objnum(3),1);
            pairnum_noncomp = testlist(testlist(:,4) == objnum(3),2);
            pairitem_noncomp = testlist(testlist(:,4) == objnum(3),3);
            pairid_noncomp = par.studytest.pairmates{category_noncomp}(pairnum_noncomp,1);
            condition_noncomp = par.studytest.pairmates{category_cue}(pairnum_noncomp,5);
            
            %array objects1
            objpos = testlist(t,7:9); %left, middle, right
            objnum_out = objnum(objpos); %for output
            
        else
            
            %for practice, just pull from the practice lists directly
            scenenum = testlist(t,1);
            objnum = testlist(t,2:4);
            
        end %end exp check
        
        %get the label for the scene
        if scenenum < 10
            snum = ['00' num2str(scenenum)];
        elseif scenenum < 100
            snum = ['0' num2str(scenenum)];
        else
            snum = num2str(scenenum);
        end
        
        if exp
            cuestim = [par.stim.labels{category_cue} '_' snum '.jpg'];
        else
            cuestim = ['prac_scene_' snum '.jpg'];
        end
        
        %get the labels for the object
        objstim = cell(1, length(objnum));
        for i = 1:length(objnum)
            
            if objnum(i) < 10
                onum = ['00' num2str(objnum(i))];
            elseif objnum(i) < 100
                onum = ['0' num2str(objnum(i))];
            else
                onum = num2str(objnum(i));
            end
            
            if exp
                objstim{i} = ['object_' onum '.jpg'];
            else
                objstim{i} = ['prac_object_' onum '.jpg'];
            end
            
        end %end i
        
        %now figure out where each of the items actually goes (left, middle, right)
        if exp
            obj1stim = objstim{objpos(1)};
            obj2stim = objstim{objpos(2)};
            obj3stim = objstim{objpos(3)};
        else
            obj1stim = objstim{1};
            obj2stim = objstim{2};
            obj3stim = objstim{3};
        end
        
        %load in the stimuli
        [cuestimvalues, ~, ~] = imread(cuestim); %cue (top)
        [obj1stimvalues, ~, ~] = imread(obj1stim);  %choice 1 (bottom left)
        [obj2stimvalues, ~, ~] = imread(obj2stim); %choice 2 (bottom middle)
        [obj3stimvalues, ~, ~] = imread(obj3stim);  %choice 3 (bottom right)
        
        %show the cue and figure out when to flip back to fixation
        Screen('PutImage', window, cuestimvalues, itemrect);
        flipped = Screen(window, 'Flip');
        ISItime = flipped + par.timing.studytest.test_cue;
        
        %flip to fixation for the ISI and figure out when to show the array
        Screen(window, 'TextSize', par.ptb.txtsize_fix);
        Screen(window, 'FillRect', par.ptb.bgcol, screenrect);
        DrawFormattedText(window,'+','center','center',par.ptb.txtcol);
        flipped = Screen(window, 'Flip',ISItime);
        arraytime = flipped + par.timing.studytest.test_ISI;
        
        %show the choice array
        Screen('PutImage', window, obj1stimvalues, obj1rect);
        Screen('PutImage', window, obj2stimvalues, obj2rect);
        Screen('PutImage', window, obj3stimvalues, obj3rect);
        Screen(window, 'Flip',arraytime);
        
        %response
        resp = [];
        respstart = GetSecs; %start of response time
        endresp = respstart + par.timing.studytest.test_array; %response deadline
        respwait = 1;
        while GetSecs < endresp
            [keyIsDown, ~, keyCode] = KbCheck; %check for key response
            if keyIsDown
                if keyCode(leftkey)
                    if isempty(resp)
                        resp = 1;
                        respstop = GetSecs;
                        respwait = 0;
                    end
                elseif keyCode(middlekey)
                    if isempty(resp)
                        resp = 2;
                        respstop = GetSecs;
                        respwait = 0;
                    end
                elseif keyCode(rightkey)
                    if isempty(resp)
                        resp = 3;
                        respstop = GetSecs;
                        respwait = 0;
                    end
                end %end of keycode check
            end  %end of kbcheck
            
            %end the trial if the participant made a response AND if they have released the key (to prevent holding down the response button)
            if respwait == 0
                [keyIsDown, ~, ~] = KbCheck;
                if ~keyIsDown
                    break;
                end
            end
            
        end %end of while
        
        %resp coding and trial RT
        if exp %only necessary for the main experiment
            if isempty(resp)
                resp = nan;
                rt = nan;
                score = nan;
            else
                rt = respstop - respstart;
                score = objpos(resp);
            end
        end %end exp check
        
        %save out the data
        if exp
            %save trial into matrix
            thistrl = [run t header.cb category_cue category_noncomp pairnum_cue pairnum_noncomp pairitem_cue pairitem_noncomp pairid_cue pairid_noncomp scenenum objnum_out condition_cue condition_noncomp objpos resp score rt];
            testdata = [testdata; thistrl]; %#ok<AGROW>
            
            %save trial to text file
            fprintf(Tid,'%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%0.3f\n',thistrl);
        end
        
        %ITI
        Screen(window, 'TextSize', par.ptb.txtsize_fix);
        Screen(window, 'FillRect', par.ptb.bgcol, screenrect);
        DrawFormattedText(window,'+','center','center',par.ptb.txtcol);
        Screen(window, 'Flip');
        WaitSecs(par.timing.studytest.test_ITI);
        
    end %end of t trials
    
    %% save your biz
    if exp
        data.studydata = studydata;
        data.testdata = testdata;
        data.endTime = fix(clock);
        save(sprintf('%s/%s_studytest_%d',header.path.subjinfo,header.subjinfo,run),'data');
    end
    
end %end of r runs

%% all done, close that mofo
Screen(window, 'TextSize', par.ptb.txtsize_inst);
DrawFormattedText(window,'End of section. Thanks!','center','center',par.ptb.txtcol);
Screen(window, 'Flip');
WaitSecs(2);
Screen('CloseAll');
Screen('Preference','SuppressAllWarnings',oldEnableFlag);
fclose('all');

end
