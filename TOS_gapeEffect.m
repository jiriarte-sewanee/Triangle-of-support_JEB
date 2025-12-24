%% TOS_gapeEffect

% This script calculates the effect of variation in muscle force
% capabilities due to gape (as a consequence of the length-tension curve)
% on JRFs and the position of the RMF with respect to the ToS.
%
% This script needs:
%   LANDMARK DATA: 
%       pp file (exported from Meshlab) with the xyz position of the
%       TMJs and the teeth
%   BONE MODEL and ROTATION DATA:
%       mat files that contains the 3D triangulated cranium and mandible
%       models and the rotation matrices for each gape angle. These data
%       are exported from the 'modelPreparation' MATLAB app
%   MUSCLE DATA:
%       mat file that contains the position, length, and moment arms data
%       for the three jaw adductor muscles at different gapes. This data is
%       obatined by running the "calculateMuscleMoments.m" script

clc;clear

%-------------------------------------------------------------------------%
% MUSCLE ARCHITECTURE (MA) DATA
n = 12 ; % number of different values for a MA variable

% Fiber length (as a fraction of muscle-tendon unit (MTU) length)
min_fiberLength = 0.5 ; % as a fraction of MTU length
max_fiberLength = 0.8 ; % as a fraction of MTU length
fiberLengths = linspace(min_fiberLength,max_fiberLength,n) ;

% Optimal length (fiber length strain from occlusion at which the muscle
% fiber produces maxmimal force)
min_optLength = 0.1 ;
max_optLength = 0.5 ;
optLength = linspace(min_optLength,max_optLength,n) ;

% pennation angle (in degrees)
pennAngle = 15 ; 

% Relative physiological cross-sectional area (PCSA)
% Assumes that all muscles produce the same relative amount of force
PCSA = [1 1 1] ;

% Create a set of a possible combinations
[FL,OL,PA] = meshgrid(fiberLengths,optLength,pennAngle) ;
tests = [FL(:) OL(:) PA(:)] ;
nTests = size(tests,1) ;

%-------------------------------------------------------------------------%
% Read the landmark data
fileInfo = dir('*.pp') ;
if isempty(fileInfo)
    error('The teeth position (.pp) file is missing')
end
tmp = readMeshlabpp(fileInfo.name) ;
points = cell2mat(tmp(:,2:end)) ;

% get the teeth position data
teethPts = points(3:end,:) ;

% get the TMJ position data
TMJpts = points(1:2,:) ;

%-------------------------------------------------------------------------%
% get the bone model data
fileInfo3 = dir('boneModels_*.mat') ;
if isempty(fileInfo3)
    error('The "boneModel" file is missing')
elseif size(fileInfo3,1)>1
    error('There is more than *one* "boneModel" file')
end
load(fileInfo3.name)

%-------------------------------------------------------------------------%
% read the muscle data
fileInfo1 = dir('MuscleMoments_*.mat') ;
if isempty(fileInfo1)
    error('The "MuscleMoments" file is missing \nPlease run the "calculateMuscleMoments.m" script')
elseif size(fileInfo1,1)>1
    error('There is more than *one* "MuscleMoments" file')
end
load(fileInfo1.name)
muscleData = muscleMoment ;

%-------------------------------------------------------------------------%
% Get the rotation data
gapeAngles = muscleData.GapeAngle ;
xTransl = muscleData.xTransl ;
zTransl = muscleData.zTransl ;
nAngles = 9 ;%numel(gapeAngles) ;

%-------------------------------------------------------------------------%
% Load the L-T curve
parentFolder = fileparts(cd) ;
load(fullfile(parentFolder,'normalized_LTcurve.mat'))

outputGape = nan(nTests,3,nAngles) ;
for q = 1:nAngles
    % create the rotation matrix for the given rotation
    rotY = troty(gapeAngles(q),'deg') ;
    T = transl(xTransl(q),0,zTransl(q)) ;
    R1 = T*rotY ;
    
    % rotate the model so that the TMJ axis is parallel to [0 1 0]
    TMJaxis2 = normalizeVector3d(diff(rotateT(TMJpts,R1))) ;
    if TMJaxis2(2)<0; TMJaxis2 = -TMJaxis2; end
    ax1 = normalizeVector3d(cross(TMJaxis2,[0 0 1])) ;
    ax3 = normalizeVector3d(cross(ax1,TMJaxis2)) ;
    R2 = eye(4) ;
    R2(1:3,1:3) = [ax1;TMJaxis2;ax3] ;
    meanTMJrot = mean(rotateT(TMJpts,R2*R1)) ;
    R2(2,4) = -meanTMJrot(2) ;
    rotT = R2*R1 ; % rotation matrix

    % rotate the teeth and TMJ points
    rotTeeth = rotateT(teethPts,rotT) ;
    rotTMJ   = rotateT(TMJpts  ,rotT) ;
    origin = mean(rotTMJ) ;

    % calculate the mean muscle resultant
    mPos = nan(3*7,3) ;
    mForce = nan(3*7,3) ;
    mForceRef = nan(3*7,3) ;

    allMusclePos = nan(nTests,3) ;
    allMuscleDir = nan(nTests,3) ;

    for w = 1:nTests % ITERATE THROUGH THE MUSCLE ARCHITECTURE CONFIGURATIONS
        % Define the muscle architecture parameters
        fLength_cf = tests(w,1) ;
        activeOffset = tests(w,2) ;
        refPennAngle = tests(w,3) ;

        for i = 1:3 % ITERATE THROUGH THE THREE MASTICATORY MUSCLES
            % get the muscle position
            musclePos = rotateT(muscleData.musclePosition{q,i},R2) ;

            % get the muscle lengths
            muscleLength = muscleData.muscleLength{q,i} ;
            muscleLengthRef = muscleData.muscleLength{1,i} ; % length at occlusion
            deltaMuscleLength = muscleLength - muscleLengthRef ;

            % calculate the fascicle length
            refFiberLength = muscleLengthRef.*fLength_cf ;

            % calculate the fiber strain
            [tmpStrain,tmpAngle] = calculateFiberStrainFromMuscleLength...
                (deltaMuscleLength,refFiberLength,refPennAngle) ;
            fiberStrain = tmpStrain.constantWidth ;
            pAngle = tmpAngle.constantWidth ;

            LTforce = interp1(LTcurve.active.Strain+activeOffset, ...
                LTcurve.active.relForce,LTcurve.active.Strain,"linear", ...
                "extrap") ;
            relForce = interp1(LTcurve.active.Strain,LTforce, ...
                fiberStrain,"linear","extrap") ;

            for j = 1:7 % ITERATE THROUGH THE MUSCLE SEGMENTS
                selMuscle = musclePos(:,3*j-2:3*j) ;
                muscleDir = -diff(selMuscle) ;

                % get the mode direction
                a = mode(muscleDir) ;
                angle = ang2vec(muscleDir(1:end-1,:),muscleDir(2:end,:),'deg') ;

                % find the last point with that direction
                selIdx = find(abs(angle)<1e-10,1,'last') ;
                if isempty(selIdx)
                    selIdx = numel(angle) ;
                end

                % determine the position and direction of the muscle force
                selMusclePos = selMuscle(selIdx+2,:) ;
                selMuscleDir = normalizeVector3d(muscleDir(selIdx,:)) ; %.*relForce(i) ;

                % create the output for muscle poisition and direction for all
                % muscles
                mPos((i-1)*7+j,:) = selMusclePos ;
                mForce((i-1)*7+j,:) = selMuscleDir.*PCSA(i).*relForce(j).*cosd(pAngle(j)) ;

                % muscle force reference. It assumes that the muscle is not
                % pennated and that its force doesn't change with gape
                mForceRef((i-1)*7+j,:) = selMuscleDir ;
            end
        end

        % calculate the position and direction of the resultant muscle
        % force by assuming that is located at Y=0 and that the lateral
        % component is canceled out
        RMF_pos = mPos.*[1 0 1] ;
        RMF = mForce.*[1 0 1] ;
        RMF_ref = mForceRef.*[1 0 1] ;

        % calculate the total force mangnitude and direction
        totalForceMag = sum(vecmag(RMF)) ;
        totalForceDir = normalizeVector3d(sum(RMF)) ;
        allMuscleDir(w,:) = totalForceDir ;

        % calculate the torque around the TMJ as the
        % dot(TMJaxis,cross(position,Force))
        d = RMF_pos-origin ;
        torques = cross(d,RMF) ;
        yTorque = torques(:,2) ;
        totalTorque = sum(yTorque) ;

        for t = 1%:size(rotTeeth,1)
            % determine the triangle of support plane
            teethPt = rotTeeth(t,:) ;
            TOSpts = triangulation([1 2 3],[rotTMJ;teethPt]) ;
            TOSplane = createPlane(TOSpts.Points) ;

            % define the balancing side joint
            if sign(rotTMJ(1,2))==sign(teethPt(2))
                pBS = rotTMJ(2,:) ;
                pWS = rotTMJ(1,:) ;
            else
                pBS = rotTMJ(1,:) ;
                pWS = rotTMJ(2,:) ;
            end

            % determine the perpendicular distance in the occlusal plane
            planeAxis = normalizeVector3d(teethPt - origin) ;
            angle1 = atan2d(planeAxis(3),planeAxis(1)) ;
            totalForce = sum(RMF) ;
            dist2Pt = totalTorque/(sind(angle1)*totalForce(1)-cosd(angle1)*totalForce(3)) ;
            ptsOnPlane = origin + dist2Pt.*[cosd(angle1) 0 sind(angle1)] ;

            % Caclculate the reaction forces at the joints
            [Fb,JRFws,JRFbs,M] = biteForce(RMF,RMF_pos,teethPt,rotTMJ) ;
            sumJRFws = sum(JRFws)./vecmag(totalForce).*100 ;
            sumJRFbs = sum(JRFbs)./vecmag(totalForce).*100 ;

            [~,JRFws_ref] = biteForce(RMF_ref,RMF_pos,teethPt,rotTMJ) ;
            totalForceRef = sum(RMF_ref) ;
            sumJRFws_ref = sum(JRFws_ref)./vecmag(totalForceRef).*100 ;

            % calculate the difference in vertical WS joint reaction force
            diffWS = sumJRFws(3) - sumJRFws_ref(3) ;

            % calculate the anterior position of the midline of the TOS by
            % estimating the intersection of the line anterior edge of the TOS
            % and the line of the sagittal midline.
            antAxisLine = createLine3d(pBS,teethPt) ;
            sagLine = createLine3d(origin,ptsOnPlane(1,:)) ;
            [~,pt1,pt2] = distanceLines3d(sagLine,antAxisLine) ;
            antTOSpt = pt1 ; % p1 should be the same as pt2

            % calculate the position of the resultant force on the occlusal
            % plane. For this, we assume that if the resultant force it's
            % at the edge of the triangle of support, then the force at the
            % ws joint will be zero, and if the muscle force is at the
            % origin, then the force should be equal and iverse to the
            % vertical force component. Then interpolate to get the point.
            sagTOSLine = createLine3d(origin,antTOSpt) ;
            resForcePos = antTOSpt + sagTOSLine(4:6).*sumJRFws(3)./100 ;
            resForcePosRef = antTOSpt + sagTOSLine(4:6).*sumJRFws_ref(3)./100 ;
            allMusclePos(w,:) = resForcePosRef ;

            % get the angle of the resultant muscle force
            invRy =inv(rotY) ;
            rotMuscleForceDir = rotateT(totalForceDir,invRy) ;
            resForceAngle = atan2d(rotMuscleForceDir(1),rotMuscleForceDir(3)) ;

            % % get the position of the anterior point of the TOS with respect to
            % % the jaw length
            % relPosAntTOS = antTOSpt(1)./jawLength ;

            % create the output file
            outputGape(w,:,q) = [resForceAngle sumJRFws(3) diffWS] ;
        end
    end
end


%-------------------------------------------------------------------------%
% PLOT THE RESULTS
%-------------------------------------------------------------------------%
set(figure(103),'windowstyle','docked','Color','w')
clf

% Determine the colormap
cm = colormap(parula(nAngles)) ;

p = panel() ;
p.pack('h',{0.2 []}) ;
p(1).pack(nAngles,1) ;
p.de.margin = 5 ;
p.margintop = 10 ;

bgColor = ones(1,3).*1 ;
xColor = 'k' ;
yColor = 'w' ;
edgeColor = ones(1,3).*0.1 ;
edgeAlpha = 0.5 ;

xLims = prctile(outputGape(:,3,:),[1 99],"all")'  ;
xLims2 = [-11.6060 11.6060] ;%[-max(abs(xLims)) max(abs(xLims))].*1.1 ;

nbins = 15 ;
for q = 1:nAngles
    % plot the distributions
    p(1,q,1).select() ;
    d1 = outputGape(:,3,q) ;
    d1(d1<xLims(1) | d1>xLims(2)) = nan ;
    histogram(d1,nbins,'FaceColor',cm(q,:), ...
        'EdgeColor','none')
    hold on
    histogram(d1,nbins,'EdgeColor',edgeColor ...
        ,'DisplayStyle','stairs') ;
    yLims = ylim ;
    axis([xLims2 yLims])
    plot([0 0],yLims,':r','LineWidth',2)
    set(gca,'Color',bgColor,'XColor',xColor,'YColor',yColor,'YTick',[])
    if q~=nAngles
        set(gca,'XTick',[])
    else
        xlabel({'Change in ws JRF','(% total muscle force)'})
    end
    text(xLims2(1)*1.1,yLims(2)*0.5, ...
        [num2str(gapeAngles(q)) ' deg'],"Color",'k', ...
        'HorizontalAlignment','right')
    if q==1
        title({'Effect of gape-force on','ws joint reaction force'},'color','k')
    end
end
p.marginbottom = 18 ;
p.margintop = 10 ;
p.fontsize = 12 ;

% save the figure
prts = strsplit(fileInfo1.name,'_') ;
expFilename = ['TOSgape_' strjoin(prts(2:4),'_') '.png'] ;
% exportgraphics(gcf,expFilename,"Resolution",300,"BackgroundColor","current")
% save('TOS_results.mat',"outputGape",'-append')