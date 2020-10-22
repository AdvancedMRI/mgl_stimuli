% adjustAnaglyphs - fix R and G guns for red-green anagyphs
%
%      usage: [ rgb ] = adjustAnaglyphs(  )
%         by: denis schluppeck
%       date: 2009-12-19
%        $Id$:
%     inputs: 
%    outputs: rgb
%
%    purpose: adjust red and green guns to make RG anaglyph goggles
%    work for L and R eye stimulation.
%
%    the first set of dots are more RED (round), 
%    the second ones more GREEN (square)
%
%    use buttons 1,2 and 3,4, respectively to adjust the nulling colours
%    DOWN and UP
%
%
function [ rgb ]=adjustAnaglyphs( varargin  )

getArgs(varargin, [], 'verbose=1');
% set some default arguments...
if ~exist('type','var'), type = 0; end
if ~exist('debug','var'), debug = 1; end
if ~exist('TR','var'), TR = 1; end
if ieNotDefined('doEyeCalib'),doEyeCalib = -1;end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other screen parameters
myscreen.autoCloseScreen = 1;
myscreen.displayname = 'projector';
myscreen.background = 'gray';
myscreen.allowpause = 1;
myscreen.saveData = 1;
myscreen.monitorGamma = 1;
myscreen.hideCursor = 1;
myscreen.TR = TR;
myscreen.eatkeys = 1;

if debug
  % gethostname and then display the stimulus on the corresponding screen
  %myscreen.screenParams{1} = {gethostname(),[],0,800,600,57,[31 23],60,1,1,1.4,[],[0 0]};
    myscreen.keyboard.nums = mglCharToKeycode({'1' '2' '3' '4'});
else
  % running at 7T for experiment
  defaultMonitorGamma = 1.8;
  myscreen.screenParams{1} = {gethostname(),'',2,1024,768,242,[69 3*69/4],60,1,1,defaultMonitorGamma,'',[0 0]}; % 7T nottingham

  % on the 7T we use buttons 9 and 8 instead of 1 and 2
  myscreen.keyboard.nums = mglCharToKeycode({'9' '8' '3' '4'  '6' '7' '0' '1' '2'});

end



% initalize the screen
myscreen = initScreen(myscreen);

% set the first task to be the fixation staircase task
% [task{1} myscreen] = fixStairInitTask(myscreen);

% set our task to have two phases. 
% one starts out with dots moving for incohrently for 10 seconds
% task{2}{1}.waitForBacktick = 1;
% task{2}{1}.seglen = 10;
% task{2}{1}.numBlocks = 1;
% task{2}{1}.parameter.dir = 0;
% task{2}{1}.parameter.coherence = 0;

% 1s "trial"
task{1}{1}.seglen = [3];
task{1}{1}.getResponse = [1];

task{1}{1}.parameter.dir = 0:60:360;
task{1}{1}.parameter.coherence = 1;
task{1}{1}.random = 1;

% initialize our task
for phaseNum = 1:length(task{1})
  [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@updateScreenCallback,@trialResponseCallback);
end

% init the stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);
stimulus = initDots(stimulus,myscreen);

% stimulus - color increment per button press
stimulus.colInc = 0.0125;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
  % update the dots
  [task{1} myscreen phaseNum] = updateTask(task{1},myscreen,phaseNum);
  % update the fixation task
  % [task{1} myscreen] = updateTask(task{1},myscreen,1);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% set output arg of adjustAnaglyphs script
rgb.dots1 = stimulus.dots.color1;
rgb.dots2 = stimulus.dots.color2;

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;
if (task.thistrial.thisseg == 1)
  stimulus.dots.coherence = task.thistrial.coherence;
else
  stimulus.dots.coherence = 0;
end
stimulus.dots.dir = task.thistrial.dir;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task, myscreen)

global stimulus
mglClearScreen;
stimulus = updateDots(stimulus,myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initDots(stimulus,myscreen)

% convert the passed in parameters to real units
if ~isfield(stimulus,'dots') || ~isfield(stimulus.dots,'rmax'), stimulus.dots.rmax = 15;,end
if ~isfield(stimulus.dots,'xcenter'), stimulus.dots.xcenter = 0;,end
if ~isfield(stimulus.dots,'ycenter'), stimulus.dots.ycenter = 0;,end
if ~isfield(stimulus.dots,'dotsize'), stimulus.dots.dotsize = 1;,end
if ~isfield(stimulus.dots,'density'), stimulus.dots.density = 3;,end
if ~isfield(stimulus.dots,'coherence'), stimulus.dots.coherence = 1;,end
if ~isfield(stimulus.dots,'speed'), stimulus.dots.speed = 5;,end
if ~isfield(stimulus.dots,'dir'), stimulus.dots.dir = 0;,end
% jwp worked the following out for his psychopy binocRiv version
% print 'red=', [1, gBalance, -1]
% print 'green', [rBalance, 0.6, -1]
% we have 0 to +1 no -1 to +1
%if ~isfield(stimulus.dots,'color1'), stimulus.dots.color1 = [1 .5  0];,end
if ~isfield(stimulus.dots,'color1'), stimulus.dots.color1 = [1  .5  .5];,end
if ~isfield(stimulus.dots,'color2'), stimulus.dots.color2 = [.5 .9 0];,end


% actually a square patch of dots that get stenciled
% so calculate width and height
stimulus.dots.width = stimulus.dots.rmax*2;
stimulus.dots.height = stimulus.dots.rmax*2;

% get the number of dots
stimulus.dots.n = round(stimulus.dots.width*stimulus.dots.height*stimulus.dots.density);

% get max and min points for dots
stimulus.dots.xmin = -stimulus.dots.width/2;
stimulus.dots.xmax = stimulus.dots.width/2;
stimulus.dots.ymin = -stimulus.dots.height/2;
stimulus.dots.ymax = stimulus.dots.height/2;

% get initial position
stimulus.dots.x = rand(1,stimulus.dots.n)*stimulus.dots.width;
stimulus.dots.y = rand(1,stimulus.dots.n)*stimulus.dots.height;

% get the step size
stimulus.dots.stepsize = stimulus.dots.speed/myscreen.framesPerSecond;

% create stencil
mglClearScreen;
mglStencilCreateBegin(1);
% and draw that oval
mglGluDisk(stimulus.dots.xcenter,stimulus.dots.ycenter,[stimulus.dots.rmax stimulus.dots.rmax]/2,[1 1 1],60);
mglStencilCreateEnd;
mglClearScreen;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to update dot positions and draw them to screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = updateDots(stimulus,myscreen)

% get the dots step
stimulus.dots.xstep = cos(pi*stimulus.dots.dir/180)*stimulus.dots.stepsize;
stimulus.dots.ystep = sin(pi*stimulus.dots.dir/180)*stimulus.dots.stepsize;

% pick a random set of dots
stimulus.dots.coherent = rand(1,stimulus.dots.n) < stimulus.dots.coherence;

% now move those dots in the right direction
stimulus.dots.x(stimulus.dots.coherent) = stimulus.dots.x(stimulus.dots.coherent)+stimulus.dots.xstep;
stimulus.dots.y(stimulus.dots.coherent) = stimulus.dots.y(stimulus.dots.coherent)+stimulus.dots.ystep;

% randomwalk rule
thisdir = rand(1,sum(~stimulus.dots.coherent))*2*pi;
stimulus.dots.x(~stimulus.dots.coherent) = stimulus.dots.x(~stimulus.dots.coherent)+cos(thisdir)*stimulus.dots.stepsize;
stimulus.dots.y(~stimulus.dots.coherent) = stimulus.dots.y(~stimulus.dots.coherent)+sin(thisdir)*stimulus.dots.stepsize;

% movshon noise
%stimulus.dots.x(~stimulus.dots.coherent) = rand(1,sum(~stimulus.dots.coherent))*stimulus.dots.width;
%stimulus.dots.y(~stimulus.dots.coherent) = rand(1,sum(~stimulus.dots.coherent))*stimulus.dots.height;

% make sure we haven't gone off the patch
% do the dots separately for left and right hand side
stimulus.dots.x(stimulus.dots.x < stimulus.dots.xmin) = stimulus.dots.x(stimulus.dots.x < stimulus.dots.xmin)+stimulus.dots.width;
stimulus.dots.x(stimulus.dots.x > stimulus.dots.xmax) = stimulus.dots.x(stimulus.dots.x > stimulus.dots.xmax)-stimulus.dots.width;
stimulus.dots.y(stimulus.dots.y < stimulus.dots.ymin) = stimulus.dots.y(stimulus.dots.y < stimulus.dots.ymin)+stimulus.dots.height;
stimulus.dots.y(stimulus.dots.y > stimulus.dots.ymax) = stimulus.dots.y(stimulus.dots.y > stimulus.dots.ymax)-stimulus.dots.height;


% draw the dots
mglStencilSelect(1);
n = stimulus.dots.n;
oneIdx = 1:floor(n/2);
twoIdx = (max(oneIdx)+1):n;

% draw half the dots in color1, the other half in color2
% used mglPoints2 before
mglGluDisk(stimulus.dots.x(oneIdx),stimulus.dots.y(oneIdx),...
	   stimulus.dots.dotsize/12,stimulus.dots.color1);

% dots
mglPoints2(stimulus.dots.x(twoIdx),stimulus.dots.y(twoIdx),...
	   stimulus.dots.dotsize*4,stimulus.dots.color2);

mglStencilSelect(0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to update colors for dots1 and dots2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [task myscreen] = trialResponseCallback(task,myscreen)

global stimulus;

% get correct or incorrect
response = find(task.thistrial.buttonState);
response = response(1);

c1 = stimulus.dots.color1;
c2 = stimulus.dots.color2;


switch response
 case 1
  %c1 = c1 + [0 -stimulus.colInc 0];
   c1 = c1 + [0 -stimulus.colInc -stimulus.colInc];
 case 2
  %c1 = c1 + [0 +stimulus.colInc 0];
   c1 = c1 + [0 +stimulus.colInc +stimulus.colInc];
 case 3
  c2 = c2 + [-stimulus.colInc 0 0];
 case 4
  c2 = c2 + [+stimulus.colInc 0 0];
end

% check that stim colors don't go beyond
c1(c1 > 1) = 1;
c1(c1 < 0) = 0;

c2(c2 > 1) = 1;
c2(c2 < 0) = 0;

stimulus.dots.color1 = c1;
stimulus.dots.color2 = c2;

fprintf('-------------------\n');
fprintf('"red dots"   c1 [%.2f %.2f %.2f]\n', c1);
fprintf('"green dots" c2 [%.2f %.2f %.2f]\n', c2)