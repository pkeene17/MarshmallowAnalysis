function data = marshmallow_nback(header,exp)
%% marshmallow_nback
% Robert Molitor, January 2020
% n-back task for for marshmallow
%
% ##input##
% header: header created from marshmallow_header
% exp: 0 = practice, 1 = main experiment
% EXAMPLE: marshmallow_prepost(header,1) will run the main task
%
% ##output##
% run
% trial
% category: 1 = beach, 2 = gazebo
% stim: stimulus number shown
% cresp: correct response
% resp: actual response
% acc: accuracy
% rt: response time (s)

%% READY, SET, GO
par = header.parameters;
data = struct('subjnum', header.subjnum);
data.header = header;
addpath(header.path.stim);
addpath(genpath('/Users/kuhllab/Desktop/Psychtoolbox'));

%% OUTPUT FILE
if exp
    outfname = sprintf('%s/%s_nback_%d',header.path.subjinfo,header.subjinfo,exp);
    if exist([outfname '.txt'],'file')
        error('The data file for this task already exists. Ending...');
    end
    
    fid=fopen([outfname '.txt'], 'w'); % open the file
    fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
        'run','trial','category','stim','cresp','resp','acc','rt'); %create text file header
end

%% SCREEN SETUP
oldEnableFlag = Screen('Preference','SuppressAllWarnings');
Screen('Preference', 'VisualDebugLevel', 1);
Screen('Preference', 'SkipSyncTests', 2); %will not show warning screen
allscreens = Screen('screens'); %how many screens are connected to the computer
screennum = max(allscreens); %choose the last one (not the main screen)
[window, screenrect] = Screen('OpenWindow',screennum,par.ptb.bgcol);
screenLeft = screenrect(RectLeft); %pixel value for left of screen
screenRight = screenrect(RectRight); %pixel value for right of screen
screenTop = screenrect(RectTop); %pixel value for top of screen
screenBottom = screenrect(RectBottom); %pixel value for bottom of screen
ycenter = (screenBottom - screenTop)/2; %ycenter pixel value
xcenter = (screenRight - screenLeft)/2; %xcenter pixel value
Screen(window, 'TextFont', 'Arial');
HideCursor;

%% IMAGE SIZING
objheight = par.stim.objheight; %height in pixels
objwidth = par.stim.objwidth; %width in pixels

%% IMAGE POSITIONING
%image displayed in the center of the screen
x1 = xcenter - objwidth/2; %left side of image
y1 = ycenter - objheight/2; %top of image
x2 = xcenter + objwidth/2; %right side of image
y2 = ycenter + objheight/2; %bottom of image
objrect = [x1 y1 x2 y2]; %rectangular coordinates for the object location

%% SET RESPONSE OPTIONS
resp1key = KbName('j'); %button to select the left object
resp2key = KbName('k'); %button to select the left object

%% RUN LOOP
if exp
    runs = size(par.nback.runs,3);
else
    runs = size(par.nback.runs_prac,3);
end

%% RUN START
data.startTime = fix(clock);

for r = 1:runs
    
    %% GET THE STIMULI
    %get the list for this run
    run = r;
    if exp %real deal
        stimlist = par.nback.runs(:,:,r);
    else %practice
        stimlist = par.nback.runs_prac(:,:,r);
    end
    
    %% START OF THE RUN
    %show instructions
    %run text first
    Screen(window, 'FillRect', par.ptb.bgcol, screenrect);
    Screen(window, 'TextSize', par.ptb.txtsize_inst);
    if exp
        starttext = sprintf('2-BACK %d/%d\n\n\n\nPress K if the current picture is the same as two pictures ago\nPress J for all other pictures\n\nPress J or K to start!',run,runs);
        
    else
        starttext = sprintf('PRACTICE 2-BACK %d/%d\n\n\n\nPress K if the current picture is the same as two pictures ago\nPress J for all other pictures\n\nPress J or K to start!',run,runs);
    end
    
    %draw the text
    DrawFormattedText(window,starttext,'center','center',par.ptb.txtcol);
    Screen(window, 'Flip');
    
    %wait for participant response to start the experiment
    respwait = 1;
    while respwait
        [keyIsDown, ~, keyCode] = KbCheck; %check for key response from the participant
        if keyIsDown
            if keyCode(resp1key) || keyCode(resp2key)
                break
            end %end of keycode check
        end  %end of kbcheck
    end %end of while
    
    Screen(window, 'TextSize', par.ptb.txtsize_fix);
    DrawFormattedText(window,'+','center','center',par.ptb.txtcol);
    Screen(window, 'Flip');
    
    %% GO THROUGH ALL THE TRIALS
    rundata = []; %trial data for analysis
    
    %start time, used to figure out stimulus onset times
    runstart = GetSecs;
    runstart_clock = fix(clock);
    
    %trial loop
    for t = 1:size(stimlist,1)
        
        %figure out when to show the next stimulus
        stimtime = runstart + par.nback.onsets(t);
        
        stimnum = stimlist(t,1);
        if stimnum < 10
            fnum = ['00' num2str(stimnum)];
        elseif stimnum < 100
            fnum = ['0' num2str(stimnum)];
        else
            fnum = num2str(stimnum);
        end
        
        if exp
            stim = [par.stim.labels{header.design.nback.category} '_' fnum '.jpg'];
        else
            stim = ['prac_scene_' fnum '.jpg'];
        end
        
        %load in the stimulus
        [objstimvalues, ~, ~] = imread(stim);
        
        %show the stimulus
        Screen('PutImage', window, objstimvalues, objrect);
        flipped = Screen(window, 'Flip', stimtime);
        
        %response
        resp = [];
        respstart = GetSecs; %start of response time
        endpic = flipped + par.timing.nback.stim;
        while GetSecs < endpic
            [keyIsDown, ~, keyCode] = KbCheck; %check for key response
            if keyIsDown
                if keyCode(resp1key)
                    if isempty(resp)
                        resp = 1;
                        respstop = GetSecs;
                    end
                elseif keyCode(resp2key)
                    if isempty(resp)
                        resp = 2;
                        respstop = GetSecs;
                    end
                end %end of keycode check
            end  %end of kbcheck
        end %end of while
        
        %remove the stimulus from the screen
        Screen(window, 'FillRect', par.ptb.bgcol, screenrect);
        Screen(window, 'TextSize', par.ptb.txtsize_fix);
        DrawFormattedText(window,'+','center','center',par.ptb.txtcol);
        Screen(window, 'Flip',endpic);
        
        %get RT and fix resp if they didn't respond in time
        if isempty(resp)
            resp = nan;
            RT = nan;
        else
            RT = respstop - respstart; %trial RT
        end
        
        %calculate accuracy
        if resp == stimlist(t,2)
            acc = 1;
        else
            acc = 0;
        end%
        
        %save trial into matrix
        thistrl = [run t header.design.nback.category stimlist(t,:) resp acc RT];
        rundata = [rundata; thistrl]; %#ok<*AGROW>
        
        %save trial to text file
        if exp
            fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%d\t%d\t%0.3f\n',thistrl);
        end
        
    end %end of t trials
    
    %wait for the last trial fixation
    WaitSecs(par.timing.nback.ITI);
    
    %accuracy by trial type
    accuracy = round(mean(rundata(:,7)) * 100);
    
    %% save your biz
    data.rundata = rundata;
    data.runstart = runstart_clock;
    data.endTime = fix(clock);
    save(sprintf('%s/%s_nback_%d_%d',header.path.subjinfo,header.subjinfo,exp,run),'data');
    
    %% BETWEEN RUN INTERVAL
    %wait a bit first
    WaitSecs(2);
    
    %display accuracy
    acctext = sprintf('Accuracy for this run: %d%%',accuracy);
    Screen(window, 'TextSize', par.ptb.txtsize_inst);
    Screen(window, 'FillRect', par.ptb.bgcol, screenrect);
    DrawFormattedText(window,acctext,'center','center',par.ptb.txtcol);
    Screen(window, 'Flip');
    WaitSecs(3);
    
    %back to fixation before the next run
    Screen(window, 'FillRect', par.ptb.bgcol, screenrect);
    Screen(window, 'TextSize', par.ptb.txtsize_fix);
    DrawFormattedText(window,'+','center','center',par.ptb.txtcol);
    Screen(window, 'Flip');
    WaitSecs(2);
    
end %end r runs

%% all done, close that mofo
Screen(window, 'TextSize', par.ptb.txtsize_inst);
DrawFormattedText(window,'End of section. Thanks!','center','center',par.ptb.txtcol);
Screen(window, 'Flip');
WaitSecs(2);
Screen('CloseAll');
Screen('Preference','SuppressAllWarnings',oldEnableFlag);
fclose('all');

end