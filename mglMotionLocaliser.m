% mglMotionLocaliser.m
%
%        $Id: mglMotionLocaliser.m,v 1.1 2016/03/14
%      usage: mglMotionLocaliser(varargin)
%         by: alex beckett
%       date: 03/14/16
%  copyright: (c) 2007 Justin Gardner (GPL see mgl/COPYING)
%    purpose: Displays a drifting dot stimulus, there are many
%             parameters you can set, you can set as many
%             or as few as you like.
%
%
%             The duty cycle
%             mglRetinotopy('dutyCycle=.15');
%             
%             THe number of cycles to run for
%             mglMotionLocaliser('numCycles=10');
%
%             The length in secs of a stimulus period
%             mglMotionLocaliser('stimulusPeriod=36');
% 
%             The number of steps that the stimulus will move in
%             one period
%             mglRetinotopy('stepsPerCycle',24');
%
%             Or instead of stimulusPeriod/stepsPerCycle one can
%             set the number of volumes per cycle and the
%             program will synch to backticks
%             mglMotionLocaliser('volumesPerCycle=24');
%
%             Eye calibration can be absent=0, at the beginning=-1
%             or at the end =1
%             mglMotionLocaliser('doEyeCalib=1');

function myscreen = mglMotionLocaliser(varargin)

% evaluate the arguments
eval(evalargs(varargin));

% setup default arguments
if exist('wedges','var') && ~isempty(wedges),wedgesOrRings = 1;end
if exist('rings','var') && ~isempty(rings),wedgesOrRings = 0;end
if ieNotDefined('wedgesOrRings'),wedgesOrRings = 1;end
if ieNotDefined('direction'),direction = 1;end
if ieNotDefined('dutyCycle'),dutyCycle = 0.5;end
if ieNotDefined('stepsPerCycle'),stepsPerCycle = 24;end
if ieNotDefined('stimulusPeriod'),stimulusPeriod = 24;end
if ieNotDefined('numCycles'),numCycles = 5;end
if ieNotDefined('doEyeCalib'),doEyeCalib = -1;end
if ieNotDefined('initialHalfCycle'),initialHalfCycle = 1;end
if ieNotDefined('debug'),debug = 0;end
if ieNotDefined('TR'),TR = 2;end

% initalize the screen
myscreen.autoCloseScreen = 1;
myscreen.allowpause = 1;
myscreen.displayname = 'projector';
% myscreen.displayname = 'onscreen'; % use this line 
myscreen.background = 'black';
myscreen.TR=3;

if debug
  % debug on eagle or laptop
  myscreen.screenParams{1} = {gethostname(),[],0,800,600,57,[31 23],60,1,1,1.4,[],[0 0]}; 
end
myscreen = initScreen(myscreen);

myscreen.keyboard.backtick = mglCharToKeycode({'5'});
% if 1%~isfield(myscreen.keyboard,'nums')
%   myscreen.keyboard.nums = mglCharToKeycode({'1' '2' '3' '4'  '6' '7' '8' '9' '0'});
myscreen.keyboard.nums = mglCharToKeycode({'9' '8' '7' '6'    '4' '3' '2' '1' '0'});
% end

% set the first task to be the fixation staircase task
global fixStimulus;
fixStimulus.fixWidth = 1;
fixStimulus.diskSize = 2;
[task{1} myscreen] = fixStairInitTask(myscreen);

% init the stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);




% set the number and length of the stimulus cycles
stimulus.numCycles = numCycles;
% set whether to have an initial half cycle
stimulus.initialHalfCycle = initialHalfCycle;
% if we have defined volumesPerCycle then get timing in volumes
if ~ieNotDefined('volumesPerCycle')
  stimulus.volumesPerCycle = volumesPerCycle;
  stimulus.stepsPerCycle = stimulus.volumesPerCycle;
% otherwise do it in seconds
else
  stimulus.stimulusPeriod = stimulusPeriod;
  % this will control how many steps the stimulus makes per cycle
  stimulus.stepsPerCycle = stepsPerCycle;
end
% set the parameters of the stimulus
% whether to display wedges or rings
% 1 for wedges/ 0 for rings
stimulus.wedgesOrRings = wedgesOrRings;
% min/max radius is the size of the stimulus
% not sure whether to hardcode the stim radius/diameter (i.e to match the
% 20 degree diameter stimuli from H+H), or whether to leave in the code
% that works out the max height/wideth and sets accordingly. Probably
% better to leave in the auto-detect, as the scanner cannot support 20
% degree stimuli (max W = 16, max H = 12)
stimulus.minRadius = 1;
stimulus.maxRadius = min(myscreen.imageWidth,myscreen.imageHeight)-1;
% stimulus.maxRadius = 20;
% direction of stimulus
stimulus.direction = direction;
% the duty cycle of the stimulus.
% For wedges, this will control the wedge size (360*dutyCycle)
% For rings, will control the ring thickness
stimulus.dutyCycle = dutyCycle;
% angle size is the size in degrees of the
% elements of the wedge that slide against each other
stimulus.elementAngleSize = 5;
% radius is the radial length of these elements
stimulus.elementRadiusSize = 2;
% radial speed of elements moving each other
stimulus.elementRadialVelocity = 7.5;
% init the stimulus

stimulus.dots.type='Opticflow';
stimulus.dots.speed=7;
stimulus.dots.theta=0;
stimulus.dots.rmax = stimulus.maxRadius;
stimulus.dots.xcenter = 0;
stimulus.dots.density = 5;
stimulus.dots.dotsize = 4;
stimulus = initDots(stimulus,myscreen);
stimulus = initWedges(stimulus,myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set our task to have a segment for each stimulus step
% and a trial for each cycle of the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
task{2}{1}.waitForBacktick = 1;
% if we are given the parametr in cycles per volume, then
% make every segment very short and let the synchToVol capture 
% the volumes to move on to the next segment 
if ~ieNotDefined('volumesPerCycle')
  task{2}{1}.seglen = 1;
  task{2}{1}.timeInVols = 1;
  % this is set so that we can end
  task{2}{1}.fudgeLastVolume = 1;
% otherwise, we are given the stimulusPeriod and stepsPerCycle and
% we compute stuff in seconds
else
%   task{2}{1}.seglen = (stimulus.stimulusPeriod/stimulus.stepsPerCycle)/n_dir_segs*ones(1,n_dir_segs);
  task{2}{1}.seglen = ones(1,stimulus.stimulusPeriod*dutyCycle);
end

task{2}{1}.parameter.coherence = 1;
task{2}{1}.parameter.speed = [1 zeros(1,(1/dutyCycle-1))];
% n_dir_segs = 4;
% task{2}{1}.dir_seglen =  (myscreen.TR)/n_dir_segs*ones(1,nsegs);

task{2}{1}.numTrials = stimulus.numCycles/dutyCycle;
% task{2}{1}.numTrials = stimulus.numCycles*stimulus.stepsPerCycle + stimulus.initialHalfCycle*round(stimulus.stepsPerCycle/2);
[task{2}{1} myscreen] = initTask(task{2}{1},myscreen,@startSegmentCallback,@updateScreenCallback,[],@startTrialCallback);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (doEyeCalib == -1)
  myscreen = eyeCalibDisp(myscreen);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task{2})) && ~myscreen.userHitEsc
  % update the dots
  [task{2} myscreen phaseNum] = updateTask(task{2},myscreen,phaseNum);
  % update the fixation task
  [task{1} myscreen] = updateTask(task{1},myscreen,1);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (doEyeCalib == 1)
  myscreen.eyecalib.prompt = 0;
  myscreen = eyeCalibDisp(myscreen);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

function [task myscreen] = startTrialCallback(task, myscreen)

global stimulus;

stimulus.dots = feval(sprintf('setDotsSpeed%s',stimulus.dots.type),stimulus.dots,stimulus.dots.speed_*task.thistrial.speed,myscreen);
%update the mask number each trial
stimulus.currentMask = stimulus.currentMask+stimulus.direction;
if stimulus.wedgesOrRings
  stimulus.currentMask = 1+mod(stimulus.currentMask-1,stimulus.wedgeN);
else
  stimulus.currentMask = 1+mod(stimulus.currentMask-1,stimulus.ringN);
end

%update the angle number each trial
stimulus.currentAngle = stimulus.currentAngle+stimulus.direction;
if stimulus.wedgesOrRings
  stimulus.currentAngle = 1+mod(stimulus.currentAngle-1,stimulus.wedgeN);
else
  stimulus.currentAngle = 1+mod(stimulus.currentAngle-1,stimulus.ringN);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;


if (task.thistrial.thisseg <= length(task.seglen))
  stimulus.coherence = task.thistrial.coherence;
  % set speed
%   stimulus.dots = feval(sprintf('setDotsSpeed%s',stimulus.dots.type),stimulus.dots,stimulus.dots.speed,myscreen);
  stimulus.dots = feval(sprintf('setDotsDir%s',stimulus.dots.type),stimulus.dots,2*mod(task.thistrial.thisseg,2)-1,myscreen);
else
  stimulus.coherence = 0;
  if strcmp(stimulus.type,'static')
    stimulus.dots = feval(sprintf('setDotsSpeed%s',stimulus.dots.type),stimulus.dots,0,myscreen);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task, myscreen)

global stimulus
mglClearScreen;
% stimulus.dots = feval(sprintf('updateDots%s',stimulus.dots.type),[stimulus.angles(stimulus.currentAngle)-stimulus.wedgeAngle stimulus.angles(stimulus.currentAngle)],stimulus.dots,myscreen);
stimulus.dots = feval(sprintf('updateDots%s',stimulus.dots.type),[0.000000001 2*pi],stimulus.dots,myscreen);
% draw the dots
if stimulus.dots.mask,mglStencilSelect(1);end
%mglPoints2(stimulus.dots.x(stimulus.dots.color==1),stimulus.dots.y(stimulus.dots.color==1),stimulus.dots.dotsize,[1 1 1]);
%mglPoints2(stimulus.dots.x(stimulus.dots.color==0),stimulus.dots.y(stimulus.dots.color==0),stimulus.dots.dotsize,[0 0 0]);
mglPoints2(stimulus.dots.x,stimulus.dots.y,stimulus.dots.dotsize,[1 1 1]);
if stimulus.dots.mask,mglStencilSelect(0);end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initWedges(stimulus,myscreen)

% round to nearest quarter of a degree, this reduces
% some edge effects
stimulus.maxRadius = floor(stimulus.maxRadius/.25)*.25;
disp(sprintf('Stimulus radius = [%0.2f %0.2f] degrees',stimulus.minRadius,stimulus.maxRadius));
% calculate some parameters
% size of wedges
stimulus.wedgeAngle = 360*stimulus.dutyCycle;
% how much to step the wedge angle by
stimulus.wedgeStepSize = 360/stimulus.stepsPerCycle;
% Based on the duty cycle calculate the ringSize.
% Note that this is not just maxRadius-minRadius * dutyCycle
% because we start the rings off the inner edge and end
% them off the outer edge (that is the screen will be blank for
% a time, rathter than showing a ring). We also have to
% compensate in our dutyCycle for this, since the effective
% duty cycle is reduced by the time that we are offscreen
% that is for two periods
dutyCycle = stimulus.dutyCycle*stimulus.stepsPerCycle/(stimulus.stepsPerCycle-2);
stimulus.ringSize = (dutyCycle*(stimulus.maxRadius-stimulus.minRadius))/(1-dutyCycle);
% get the radii for the rings
minRadius = stimulus.minRadius-stimulus.ringSize;
maxRadius = stimulus.maxRadius+stimulus.ringSize;
% now we have the inner and outer ring radius that will be used
% add a little fudge factor so that we don't get any rings
% with only a small bit showing
epsilon = 0.1;
stimulus.ringRadiusMin = max(0,minRadius:(epsilon+stimulus.maxRadius-minRadius)/(stimulus.stepsPerCycle-1):stimulus.maxRadius+epsilon);
stimulus.ringRadiusMax = stimulus.minRadius:(maxRadius-stimulus.minRadius)/(stimulus.stepsPerCycle-1):maxRadius;

% we only need to recompute the mglQuad points of the elements if something has
% changed in the stimulus
if ~isfield(stimulus,'last') || ~isfield(stimulus,'x') || ...
  (stimulus.elementAngleSize ~= stimulus.last.elementAngleSize) || ...
  (stimulus.elementRadiusSize ~= stimulus.last.elementRadiusSize) || ...
  (stimulus.elementRadialVelocity ~= stimulus.last.elementRadialVelocity) || ...
  (stimulus.maxRadius ~= stimulus.last.maxRadius) || ...
  (stimulus.minRadius ~= stimulus.last.minRadius)
  % all the angles that the elements will be made up of
  allAngles = (0:stimulus.elementAngleSize:(360-stimulus.elementAngleSize));
  % all the phases. The phase refers to the radial position of the
  % black and white pattern (the pattern that is seen as moving
  % in the stimulus). There are two sets here since the wedges slide
  % against each other. That is every other sector will go in a 
  % different direction. 
  allPhases1 = 0:(stimulus.elementRadialVelocity/myscreen.framesPerSecond):(stimulus.elementRadiusSize*2);
  allPhases2 = fliplr(allPhases1);
  disppercent(-inf,'(mglRetinotopy) Calculating coordinates of elements in stimulus pattern');
  for phaseNum = 1:length(allPhases1)
    stimulus.x{phaseNum} = [];stimulus.y{phaseNum} = [];stimulus.c{phaseNum} = [];
    for angleNum = 1:length(allAngles)
      % get the angle
      angle = allAngles(angleNum);
      % choose which phase we are going to be
      if isodd(angleNum)
	thisMinRadius = stimulus.minRadius-allPhases1(phaseNum);
      else
	thisMinRadius = stimulus.minRadius-allPhases2(phaseNum);
      end
      % all the radiuses
      allRadius = thisMinRadius:stimulus.elementRadiusSize:stimulus.maxRadius;
      % now create all the quads for this wedge
      for radiusNum = 1:length(allRadius)
	radius = allRadius(radiusNum);
	if (radius+stimulus.elementRadiusSize) >= stimulus.minRadius
	  radius1 = max(radius,stimulus.minRadius);
	  radius2 = min(radius+stimulus.elementRadiusSize,stimulus.maxRadius);
	  % calculate in polar angle coordinates the corners of this quad
	  r = [radius1 radius1 radius2 radius2];
	  a = [angle angle+stimulus.elementAngleSize angle+stimulus.elementAngleSize angle];
	  % convert into rectilinear coordinates and save in array
	  stimulus.x{phaseNum}(:,end+1) = r.*cos(d2r(a));
	  stimulus.y{phaseNum}(:,end+1) = r.*sin(d2r(a));
	  % also calculate what color we ant
	  stimulus.c{phaseNum}(:,end+1) = [1 1 1]*(isodd(radiusNum+isodd(angleNum)));
	end
      end
    end
    disppercent(phaseNum/length(allPhases1));
  end
  disppercent(inf);
  stimulus.n = length(allPhases1);
  stimulus.phaseNum = 1;
else
  disp(sprintf('(mglRetinotopy) Using precomputed stimulus pattern'));
end
% remember these parameters, so that we can know whether we
% need to recompute
stimulus.last.elementRadiusSize = stimulus.elementRadiusSize;
stimulus.last.elementAngleSize = stimulus.elementAngleSize;
stimulus.last.elementRadialVelocity = stimulus.elementRadialVelocity;
stimulus.last.maxRadius = stimulus.maxRadius;
stimulus.last.minRadius = stimulus.minRadius;

% new we calculate the masks that cover the stimulus so that we can
% have either rings or wedges, we start by making a set of wedge masks
angles = (0:stimulus.wedgeStepSize:(360-stimulus.wedgeStepSize))+90+stimulus.wedgeAngle/2;
stimulus.angles = angles*pi/180;
stimulus.angles(stimulus.angles>2*pi) = stimulus.angles(stimulus.angles>2*pi)-2*pi;
% create masks for wedges
for angleNum = 1:length(angles)
  angle = angles(angleNum);
  % init the wedge mask values
  stimulus.maskWedgeX{angleNum} = [];
  stimulus.maskWedgeY{angleNum} = [];
  % create a polygon that spares the wedge that we want
  % start in the center, compute it in radial coordinates
  r = 0;a = 0;
  % and go around the angles except for the wedge we want
  for vertexAngle = angle:(angle+360-stimulus.wedgeAngle);
    r(end+1) = stimulus.maxRadius+1;
    a(end+1) = vertexAngle;
  end
  % and end up in the center
  r(end+1) = 0;
  a(end+1) = 0;
  % now convert to rectilinear
  stimulus.maskWedgeX{angleNum}(:,end+1) = r.*cos(d2r(a));
  stimulus.maskWedgeY{angleNum}(:,end+1) = r.*sin(d2r(a));
end
stimulus.wedgeN = length(angles);


stimulus.wedgeAngle = stimulus.wedgeAngle*pi/180;
% now we will make the masks for the rings. We will
% make an inner and outer set of ring masks so that we
% can cover up everything but one ring of the stimulus
for radiusNum = 1:length(stimulus.ringRadiusMin)
  % create the inner mask
  stimulus.maskInnerX{radiusNum} = [];
  stimulus.maskInnerY{radiusNum} = [];
  % compute in radial coordinates
  r = [0];a = [0];
  for angle = 0:stimulus.elementAngleSize:360
    r(end+1) = stimulus.ringRadiusMin(radiusNum);
    a(end+1) = angle;
  end
  r(end+1) = 0;a(end+1) = 0;
  % now convert to rectilinear
  stimulus.maskInnerX{radiusNum}(:,end+1) = r.*cos(d2r(a));
  stimulus.maskInnerY{radiusNum}(:,end+1) = r.*sin(d2r(a));
  % create the outer mask, this will be 
  % a set of quads that make a torus
  stimulus.maskOuterX{radiusNum} = [];
  stimulus.maskOuterY{radiusNum} = [];
  allAngles = 0:stimulus.elementAngleSize:360;
  for angleNum = 1:length(allAngles)
    angle = allAngles(angleNum);
    r = stimulus.ringRadiusMax(radiusNum);
    a = angle;
    r(end+1) = stimulus.maxRadius+1;
    a(end+1) = angle;
    r(end+1) = stimulus.maxRadius+1;
    a(end+1) = angle+stimulus.elementAngleSize;
    r(end+1) = stimulus.ringRadiusMax(radiusNum);
    a(end+1) = angle+stimulus.elementAngleSize;
    % convert to rectilinear
    stimulus.maskOuterX{radiusNum}(:,angleNum) = r.*cos(d2r(a));
    stimulus.maskOuterY{radiusNum}(:,angleNum) = r.*sin(d2r(a));
    stimulus.maskOuterC{radiusNum}(:,angleNum) = [0.5 0.5 0.5];
  end
end
stimulus.ringN = length(stimulus.ringRadiusMin);
  
% set the current mask that will be displayed
if stimulus.direction == 1
  stimulus.currentMask = 0;
else
  stimulus.currentMask = 2;
end

% if we are supposed to start halfway through
if stimulus.initialHalfCycle
  stimulus.currentMask = stimulus.currentMask+round(stimulus.stepsPerCycle/2);
end


% set the start angle that will be displayed
if stimulus.direction == 1
  stimulus.currentAngle = 0;
else
  stimulus.currentAngle = 2;
end

% if we are supposed to start halfway through
if stimulus.initialHalfCycle
  stimulus.currentAngle = stimulus.currentAngle+round(stimulus.stepsPerCycle/2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to draw wedges to screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = updateWedges(stimulus,myscreen)

% update the phase of the sliding wedges
stimulus.phaseNum = 1+mod(stimulus.phaseNum,stimulus.n);
% draw the whole stimulus pattern
mglQuad(stimulus.x{stimulus.phaseNum},stimulus.y{stimulus.phaseNum},stimulus.c{stimulus.phaseNum},1);

% mask out to get a wedge
if stimulus.wedgesOrRings
  mglPolygon(stimulus.maskWedgeX{stimulus.currentMask},stimulus.maskWedgeY{stimulus.currentMask},[0.5 0.5 0.5]);
% or mask out to get a ring
else
  mglPolygon(stimulus.maskInnerX{stimulus.currentMask},stimulus.maskInnerY{stimulus.currentMask},[0.5 0.5 0.5]);
  mglQuad(stimulus.maskOuterX{stimulus.currentMask},stimulus.maskOuterY{stimulus.currentMask},stimulus.maskOuterC{stimulus.currentMask});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initDots(stimulus,myscreen)


% convert the passed in parameters to real units
if ~isfield(stimulus,'dots') || ~isfield(stimulus.dots,'rmax'), stimulus.dots.rmax = min(myscreen.imageWidth,myscreen.imageHeight);,end
if ~isfield(stimulus.dots,'xcenter'), stimulus.dots.xcenter = 0;,end
if ~isfield(stimulus.dots,'ycenter'), stimulus.dots.ycenter = 0;,end
if ~isfield(stimulus.dots,'dotsize'), stimulus.dots.dotsize = 4;,end
if ~isfield(stimulus.dots,'density'), stimulus.dots.density = 5;,end
if ~isfield(stimulus.dots,'coherence'), stimulus.dots.coherence = 0;,end
if ~isfield(stimulus.dots,'speed'), stimulus.dots.speed = stimulus.speed;,end
if ~isfield(stimulus.dots,'dir'), stimulus.dots.dir = 0;,end
if ~isfield(stimulus.dots,'angle'), stimulus.dots.angle = 90;,end
stimulus.dots.mask = 1;
if ~isfield(stimulus.dots,'mirror'),stimulus.dots.mirror.on = 0;end
stimulus.dots.speed_=stimulus.dots.speed;
% stimulus.dots.rmax=25;


% update the dots
stimulus.dots = feval(sprintf('initDots%s',stimulus.dots.type),stimulus.dots,myscreen);

% set color
stimulus.dots.color = ones(stimulus.dots.n,1);
%stimulus.dots.color(rand(1,stimulus.dots.n)>0.5) = 1;

% create stencil
if stimulus.dots.mask
  mglClearScreen;
  mglStencilCreateBegin(1);
  % and draw that oval
  mglGluDisk(stimulus.dots.xcenter,stimulus.dots.ycenter,[stimulus.dots.rmax stimulus.dots.rmax]/2,[1 1 1],60);
  mglStencilCreateEnd;
  mglClearScreen;
end



%----------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%stuff for opticflow dots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots for opticflow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updateDotsOpticflow(angles,dots,myscreen)
 
[wedge_x wedge_y] = pol2cart(angles,[100 100]);

wedge_x = wedge_x/myscreen.imageWidth;
wedge_y = wedge_y/myscreen.imageHeight;

proj_angles=atan2(wedge_y,wedge_x);
proj_angles(proj_angles<0) = proj_angles(proj_angles<0)+2*pi;

phi = atan2(dots.yproj, dots.xproj);
phi(phi<0) = phi(phi<0)+2*pi;

r = sqrt(dots.yproj.^2+dots.xproj.^2);
if proj_angles(1) > proj_angles(2)
   idx = (phi >= proj_angles(1)) | (phi <= proj_angles(2));
else
    idx = (phi >= proj_angles(1)) & (phi <= proj_angles(2));
end


% get the coherent and incoherent dots
%if (dots.coherency ~= coherence)
coherence=1;
  dots.incoherent = rand(1,dots.n) > coherence;
  dots.incoherentn = sum(dots.incoherent);
  dots.coherent = ~dots.incoherent;
  dots.coherency = coherence;
  % generate a random transformation matrix for each incoherent point
  dots.randT = rand(3,dots.incoherentn)-0.5;
  % and normalize the transformation to have the same length
  % (i.e. speed) as the real transformation matrix
  dots.randT = sqrt(sum(dots.T.^2))*dots.randT./([1 1 1]'*sqrt(sum(dots.randT.^2)));
%end

% update relative position of dots in 3-space to observer
dots.X(dots.coherent) = dots.X(dots.coherent)-dots.T(1);
dots.Y(dots.coherent) = dots.Y(dots.coherent)-dots.T(2);
dots.Z(dots.coherent) = dots.Z(dots.coherent)-dots.T(3);
dots.Z_=dots.Z_-dots.T(3);

% now move the incoherent points according to the random trasnformation
dots.X(dots.incoherent) = dots.X(dots.incoherent)-dots.randT(1,:);
dots.Y(dots.incoherent) = dots.Y(dots.incoherent)-dots.randT(2,:);
dots.Z(dots.incoherent) = dots.Z(dots.incoherent)-dots.randT(3,:);

% get all points that have fallen off the screen
offscreen = dots.Z<dots.minZ;

% and put them at the furthest distance
dots.Z(offscreen) = dots.Z(offscreen)-dots.minZ+dots.maxZ;

% get all points that have fallen out of view
offscreen = dots.Z>dots.maxZ;

% and move them to the front plane

dots.Z(offscreen) =dots.Z(offscreen)-dots.maxZ+dots.minZ;

% put points fallen off the X edge back
offscreen = dots.X < -dots.maxX;
dots.X(offscreen) = dots.X(offscreen)+2*dots.maxX;
offscreen = dots.X > dots.maxX;
dots.X(offscreen) = dots.X(offscreen)-2*dots.maxX;

% put points fallen off the Y edge back
offscreen = dots.Y < -dots.maxY;
dots.Y(offscreen) = dots.Y(offscreen)+2*dots.maxY;
offscreen = dots.Y > dots.maxY;
dots.Y(offscreen) = dots.Y(offscreen)-2*dots.maxY;



% project on to screen
dots.xproj = dots.f*dots.X./dots.Z;
dots.yproj = dots.f*dots.Y./dots.Z;

% restrict updating to dots that are currently within the MOVING wedge
% make phi between 0 and 2pi

if sum(idx)<(dots.n/24),keyboard,end
dots.x(idx) = dots.xproj(idx)*myscreen.imageWidth;
dots.y(idx) = dots.yproj(idx)*myscreen.imageHeight;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the dots speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = setDotsSpeedOpticflow(dots,speed,myscreen)

% get the step size
dots.speed = speed;
dots.T = [0 0 dots.speed/myscreen.framesPerSecond];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the dots direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = setDotsDirOpticflow(dots,direction,myscreen)

% get the step size
dots.T = [0 0 direction*dots.speed/myscreen.framesPerSecond];
dots.R = makerotmatrix(dots.theta*direction);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create dots for optic flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initDotsOpticflow(dots,myscreen)

% focal length to projection plane
% projection plane is defined to be 
% 1 unit wide and high, so with 
% this focal length, we are looking at
% a view of the world with a 90 deg fov
dots.f = .5;

% translation and rotation matrices
dots.T = [0 0 dots.speed/myscreen.framesPerSecond];
dots.R = makerotmatrix(dots.theta);

% maximum depth of points
dots.maxZ = 10;dots.minZ = dots.f;
dots.maxX = 10;
dots.maxY = 10;

% make a brick of points
dots.n = round(myscreen.imageWidth*myscreen.imageHeight*dots.density);

% initial position of dots
dots.X = 2*dots.maxX*rand(1,dots.n)-dots.maxX;
dots.Y = 2*dots.maxY*rand(1,dots.n)-dots.maxY;
dots.Z = (dots.maxZ-dots.minZ)*rand(1,dots.n)+dots.minZ;
dots.Z_=dots.Z(1);
% dots.Z = 8*ones(1,dots.n)+dots.minZ;

% get projection on to plane
dots.xproj = dots.f*dots.X./dots.Z;
dots.yproj = dots.f*dots.Y./dots.Z;

% put into screen coordinates
dots.x = dots.xproj*myscreen.imageWidth;
dots.y = dots.yproj*myscreen.imageHeight;

% set incoherent dots to 0
dots.coherency = 1;
dots.incoherent = rand(1,dots.n) > dots.coherency;
dots.incoherentn = sum(dots.incoherent);
dots.coherent = ~dots.incoherent;

dots.randT = zeros(3,dots.incoherentn);






%----------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%stuff for linear dots (can't use without compresion issues)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the dots speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = setDotsSpeedLinear(dots,speed,myscreen)

% get the step size
dots.speed = speed;
dots.T = [0 0 dots.speed/myscreen.framesPerSecond];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the dots direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = setDotsDirLinear(dots,direction,myscreen)

% get the step size
if direction==1
    dots.Z=dots.Z_;
end
dots.T = [0 0 direction*dots.speed/myscreen.framesPerSecond];
dots.R = makerotmatrix(dots.theta*direction);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots for linear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updateDotsLinear(angles,dots,myscreen)

[wedge_x wedge_y] = pol2cart(angles,[100 100]);

wedge_x = wedge_x/myscreen.imageWidth;
wedge_y = wedge_y/myscreen.imageHeight;

proj_angles=atan2(wedge_y,wedge_x);
proj_angles(proj_angles<0) = proj_angles(proj_angles<0)+2*pi;

phi = atan2(dots.yproj, dots.xproj);
phi(phi<0) = phi(phi<0)+2*pi;

r = sqrt(dots.yproj.^2+dots.xproj.^2);
if proj_angles(1) > proj_angles(2)
   idx = (phi >= proj_angles(1)) | (phi <= proj_angles(2));
else
    idx = (phi >= proj_angles(1)) & (phi <= proj_angles(2));
end


% get the coherent and incoherent dots
%if (dots.coherency ~= coherence)
coherence=1;
  dots.incoherent = rand(1,dots.n) > coherence;
  dots.incoherentn = sum(dots.incoherent);
  dots.coherent = ~dots.incoherent;
  dots.coherency = coherence;
  % generate a random transformation matrix for each incoherent point
  dots.randT = rand(3,dots.incoherentn)-0.5;
  % and normalize the transformation to have the same length
  % (i.e. speed) as the real transformation matrix
  dots.randT = sqrt(sum(dots.T.^2))*dots.randT./([1 1 1]'*sqrt(sum(dots.randT.^2)));
%end

% update relative position of dots in 3-space to observer
dots.X(dots.coherent) = dots.X(dots.coherent)-dots.T(1);
dots.Y(dots.coherent) = dots.Y(dots.coherent)-dots.T(2);
dots.Z(dots.coherent) = dots.Z(dots.coherent)-dots.T(3);

% now move the incoherent points according to the random trasnformation
dots.X(dots.incoherent) = dots.X(dots.incoherent)-dots.randT(1,:);
dots.Y(dots.incoherent) = dots.Y(dots.incoherent)-dots.randT(2,:);
dots.Z(dots.incoherent) = dots.Z(dots.incoherent)-dots.randT(3,:);

% get all points that have fallen off the screen
offscreen = dots.Z<dots.minZ;

% and put them at the furthest distance
dots.Z(offscreen) = dots.Z(offscreen)-dots.minZ+dots.maxZ;

% get all points that have fallen out of view
offscreen = dots.Z>dots.maxZ;

% and move them to the front plane

dots.Z(offscreen) =dots.Z(offscreen)-dots.maxZ+dots.minZ;

% put points fallen off the X edge back
offscreen = dots.X < -dots.maxX;
dots.X(offscreen) = dots.X(offscreen)+2*dots.maxX;
offscreen = dots.X > dots.maxX;
dots.X(offscreen) = dots.X(offscreen)-2*dots.maxX;

% put points fallen off the Y edge back
offscreen = dots.Y < -dots.maxY;
dots.Y(offscreen) = dots.Y(offscreen)+2*dots.maxY;
offscreen = dots.Y > dots.maxY;
dots.Y(offscreen) = dots.Y(offscreen)-2*dots.maxY;



% project on to screen
dots.xproj = dots.f*dots.X./dots.Z;
dots.yproj = dots.f*dots.Y./dots.Z;

% restrict updating to dots that are currently within the MOVING wedge
% make phi between 0 and 2pi

if sum(idx)<(dots.n/24),keyboard,end
dots.x(idx) = dots.xproj(idx)*myscreen.imageWidth;
dots.y(idx) = dots.yproj(idx)*myscreen.imageHeight;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create dots for linear2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initDotsLinear(dots,myscreen)

% focal length to projection plane
% projection plane is defined to be 
% 1 unit wide and high, so with 
% this focal length, we are looking at
% a view of the world with a 90 deg fov
dots.f = .5;

% translation and rotation matrices
dots.T = [0 0 dots.speed/myscreen.framesPerSecond];
dots.R = makerotmatrix(dots.theta);

% maximum depth of points
dots.maxZ = 10;dots.minZ = dots.f;
dots.maxX = 10;
dots.maxY = 10;

% make a brick of points
dots.n = round(myscreen.imageWidth*myscreen.imageHeight*dots.density);

% initial position of dots
dots.X = 2*dots.maxX*rand(1,dots.n)-dots.maxX;
dots.Y = 2*dots.maxY*rand(1,dots.n)-dots.maxY;
% dots.Z = (dots.maxZ-dots.minZ)*rand(1,dots.n)+dots.minZ;
dots.Z = 8*ones(1,dots.n);
dots.Z_=dots.Z;

% get projection on to plane
dots.xproj = dots.f*dots.X./dots.Z;
dots.yproj = dots.f*dots.Y./dots.Z;

% put into screen coordinates
dots.x = dots.xproj*myscreen.imageWidth;
dots.y = dots.yproj*myscreen.imageHeight;

% set incoherent dots to 0
dots.coherency = 1;
dots.incoherent = rand(1,dots.n) > dots.coherency;
dots.incoherentn = sum(dots.incoherent);
dots.coherent = ~dots.incoherent;

dots.randT = zeros(3,dots.incoherentn);

%----------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to calculate rotation, not being used currently
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m = makerotmatrix(theta)

% check command line arguments
if (nargin ~= 1)
  help rotmatrix;
  return
end

theta = d2r(theta);
m = [cos(theta) sin(theta);-sin(theta) cos(theta)];

