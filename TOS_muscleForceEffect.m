%% TOS_muscleForceEffect

% This scripts calculates the effect of differences in relative force of
% the three jaw adductor muscles on joint reaction forces and the position
% of the resultant muscle force with respect to the triangle of support
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

clear; clc

%-------------------------------------------------------------------------%
% RELATIVE MUSCLE FORCES 
%-------------------------------------------------------------------------%
% Define a series of relative muscle forces for the masseter, medial
% pterygoid, and temporalis. The sum for each set of relative muscle forces
% is 1.

% create a grid of relative muscle forces for the masseter and for the
% medial pterygoid
mForce_min = 0.2 ;
mForce_max = 0.6 ;
nMuscleForces = 17 ;
MPt = linspace(mForce_min,mForce_max,nMuscleForces)' ;
Mass  = linspace(mForce_min,mForce_max,nMuscleForces)' ;
[mMass,mMPte] = meshgrid(Mass,MPt) ;

muscleForces = [mMass(:) mMPte(:)] ;
muscleForces(:,3) = 1-sum(muscleForces,2) ;
muscleForces(muscleForces(:,3)<mForce_min-eps,:) = [] ;
muscleForces(end+1,:) = [1/3 1/3 1/3] ;
nMuscleForceSets = size(muscleForces,1) ;

% find the indices of the muscle forces that will be used as reference
% equal muscle forces
idx_equalMuscles = nMuscleForceSets ;

% max Masseter
tmpMuscles = muscleForces(muscleForces(:,1)==mForce_max,:) ;
[~,idx] = min(mean(tmpMuscles(:,2:3))) ;
idx_maxMass = find(all(muscleForces==tmpMuscles(idx,:),2)) ;

% max Medial Pterygoid
tmpMuscles = muscleForces(muscleForces(:,2)==mForce_max,:) ;
[~,idx] = min(mean(tmpMuscles(:,[1 3]))) ;
idx_maxMPte = find(all(muscleForces==tmpMuscles(idx,:),2)) ;

% max Medial Pterygoid
tmpMuscles = muscleForces(muscleForces(:,3)>mForce_max-eps,:) ;
[~,idx] = min(mean(tmpMuscles(:,1:2))) ;
idx_maxTemp = find(all(muscleForces==tmpMuscles(idx,:),2)) ;


%-------------------------------------------------------------------------%
% READ THE TMJ AND TEETH POSITION DATA 
%-------------------------------------------------------------------------%
% position data comes from data exported from Meshlab as .pp file

% Read the points data
fileInfo = dir('*.pp') ;
if isempty(fileInfo)
    error('The teeth position (.pp) file is missing')
end
tmp = readMeshlabpp(fileInfo.name) ;

points = cell2mat(tmp(:,2:end)) ;

% get the teeth data (ordered from posterior to anterior)
teethPts = points(3:end,:) ;

% get the TMJ data
TMJpts = points(1:2,:) ;
meanTMJpts = mean(TMJpts) ;

%-------------------------------------------------------------------------%
% READ BONE MODEL AND MUSCLE POSITION DATA 
%-------------------------------------------------------------------------%

% get the bone model data
fileInfo3 = dir('boneModels_*.mat') ;
if isempty(fileInfo3)
    error('The "boneModel" file is missing')
elseif size(fileInfo3,1)>1
    error('There is more than *one* "boneModel" file')
end
load(fileInfo3.name)

% read the muscle data
fileInfo1 = dir('MuscleMoments_*.mat') ;
if isempty(fileInfo1)
    error('The "MuscleMoments" file is missing \nPlease run the "calculateMuscleMoments.m" script')
elseif size(fileInfo1,1)>1
    error('There is more than *one* "MuscleMoments" file')
end
load(fileInfo1.name)

%-------------------------------------------------------------------------%
% READ MODEL ROTATION DATA
%-------------------------------------------------------------------------%

% read the rotations data
rotations = muscleMoment(:,1:3) ;
gapeAngles = rotations.GapeAngle ;
nAngles = numel(gapeAngles) ; 
R1 = transl(rotations{1,2},0,rotations{1,3})*troty(rotations{1,1}) ;

%-------------------------------------------------------------------------%
% PREPARE THE MODEL FOR BIOMECHANICAL CALCULATIONS
%-------------------------------------------------------------------------%

% get the points of the occlusal plane. This will be the position of the
% last molar and the TMJs projected at the last molar Z value.
occlPts = [TMJpts;teethPts(1,:)] ;
occlPts(1:2,3) = teethPts(1,3) + (TMJpts(:,3)-meanTMJpts(3)) ;

%-------------------------------------------------------------------------%
% Initialize figure
%-------------------------------------------------------------------------%
set(figure(101),'windowstyle','docked','Color','w')
clf

cm = colormap(parula(nAngles)) ;
global ColorOrder, ColorOrder = [] ;

%-------------------------------------------------------------------------%
% CALCULATE BIOMECHANICAL PARAMETERS
%-------------------------------------------------------------------------%

% initialize the output files
outputMuscles = nan(nMuscleForceSets,9,nAngles) ;
muscleEffects = nan(nAngles,3) ;

% Iterate through the different gape angles
for q = 1:nAngles

    % create the rotation matrix for the given rotation
    rotY = troty(gapeAngles(q),'deg') ; 
    T = transl(rotations.xTransl(q),0,rotations.zTransl(q)) ; 
    R1 = T*rotY ; % rotation matrix
    
    % rotate the model so that the TMJ axis (the axis that connects both
    % TMJs) is aligned with the Y-axis
    TMJaxis2 = normalizeVector3d(diff(rotateT(TMJpts,R1))) ;
    if TMJaxis2(2)<0; TMJaxis2 = -TMJaxis2; end
    ax1 = normalizeVector3d(cross(TMJaxis2,[0 0 1])) ;
    ax3 = normalizeVector3d(cross(ax1,TMJaxis2)) ;
    R2 = eye(4) ;
    R2(1:3,1:3) = [ax1;TMJaxis2;ax3] ; % R2 rotates the already rotated values so that the TMJ axis is parallel to the horizontal
    meanTMJrot = rotateT(meanTMJpts,R2*R1) ;
    R2(2,4) = -meanTMJrot(2) ;
    RM = R2*R1 ; % final rotation matrix

    % rotate the teeth and TMJ points
    rotTeeth = rotateT(teethPts,RM) ;
    rotTMJ   = rotateT(TMJpts  ,RM) ;
    origin = mean(rotTMJ) ;
    rotOcclPts = rotateT(occlPts,RM) ;

    % % calculate the mean muscle resultant

    % Initilize the muscle force data
    mPos = nan(3*7,3) ; % muscle position
    mForce = nan(3*7,3) ; % muscle force vector
    mPosSingle = nan(3,3) ; 
    mForceSingle = nan(3,3) ;

    % Initialize the output of the resultant muscle forces
    allMusclePosOnToS = nan(nMuscleForceSets,3) ;
    allMusclePosOnOccl = nan(nMuscleForceSets,3) ;
    allMuscleDir = nan(nMuscleForceSets,3) ;

    % get the muscle data for the given gape angle -----------------------%
    data = muscleMoment(q,:) ;

    % Iterate through the different sets of relative muscle forces
    for w = 1:nMuscleForceSets

        % define the relative muscle force to test
        relForce = muscleForces(w,:) ;

        % Iterate though the three jaw adductor muscles
        nMuscles = 3 ;
        for i = 1:nMuscles

            % get the muscle position for the given angle
            tmpMusclePos = data.musclePosition{i} ;

            % rotate the muscle position to match the model with aligned
            % TMJ axis
            musclePos = rotateT(tmpMusclePos,R2) ;

            % ITERATE THROUGH MUSCLE SEGMENTS
            nSegments = 7 ;
            for j = 1:nSegments

                % get the position of the selected muscle segment
                selMusclePos = musclePos(:,3*j-2:3*j) ;

                % get the mode direction
                muscleDir = -diff(selMusclePos) ;
                angleSegment = ang2vec(muscleDir(1:end-1,:),muscleDir(2:end,:),'deg') ;

                % find the last point with that direction
                selIdx = find(abs(angleSegment)<1e-10,1,'last') ;
                if isempty(selIdx)
                    selIdx = numel(angleSegment) ;
                end

                % determine the position and direction of the muscle force
                selPos = selMusclePos(selIdx+2,:) ;
                selForce = normalizeVector3d(muscleDir(selIdx,:)).*relForce(i) ;

                % create the output muscle force if we want to model the
                % force each muscle as one single muscle segment
                if j==4 % if it is the middle segment
                    mPosSingle(i,:) = selPos ;
                    mForceSingle(i,:) = selForce ;
                end

                % create the output for muscle position and direction for all
                % muscles
                mPos((i-1)*7+j,:) = selPos ;
                mForce((i-1)*7+j,:) = selForce ;
            end
        end

        % because the muscle position was determined for just one side of
        % the mandible, we assume that the resultant muscle force from both
        % left and right sides muscles will have the same direction but
        % with their Y component canceled out
        mResultantPos = mPos.*[1 0 1] ;
        mResultantForce = mForce.*[1 0 1] ;
        mResultantPos_single = mPosSingle.*[1 0 1] ;
        mResultantForce_single = mForceSingle.*[1 0 1] ;

        % calculate the total force mangnitude and direction
        totalForceDir = normalizeVector3d(sum(mResultantForce)) ;
        allMuscleDir(w,:) = totalForceDir ;

        % get the angle of the resultant muscle force
        resForceAngle = atan2d(totalForceDir(1),totalForceDir(3)) - gapeAngles(q) ;

        % calculate the torque around the origin as the
        % cross(position,Force)
        d = mResultantPos - origin ;
        torques = cross(d,mResultantForce) ;
        yTorque = torques(:,2) ;
        totalTorque = sum(yTorque) ;

        % ITERATE THROUGH THE DIFFERENT BITING POINTS
        % but for now we would just calculate the forces at the most
        % posterior teeth
        for t = 1%:size(rotTeeth,1)

            % determine the ToS plane
            teethPt = rotTeeth(t,:) ;
            TOSpts = triangulation([1 2 3],[rotTMJ;teethPt]) ;
            occlPlane = triangulation([1 2 3],rotOcclPts) ;

            % define the balancing side joint
            if sign(rotTMJ(1,2))==sign(teethPt(2))
                pBS = rotTMJ(2,:) ;
                pWS = rotTMJ(1,:) ;
            else
                pBS = rotTMJ(1,:) ;
                pWS = rotTMJ(2,:) ;
            end

            % calculate the position of the force vector on the ToS plane
            % that produces the same moments and reaction force. This would
            % be where the resultant force is, but it doesn't simply
            % explain the mechanics of the system.
            planeAxis = normalizeVector3d(teethPt - origin) ;
            ToSangle = atan2d(planeAxis(3),planeAxis(1)) ;
            totalForce = sum(mResultantForce) ;
            dist2Pt = totalTorque/(sind(ToSangle)*totalForce(1)-cosd(ToSangle)*totalForce(3)) ;
            ptOnToSplane = origin + dist2Pt.*[cosd(ToSangle) 0 sind(ToSangle)] ;
            allMusclePosOnToS(w,:) = ptOnToSplane ;

            % calculate the anterior position of the midline of the TOS by
            % estimating the intersection of the line anterior edge of the TOS
            % and the line of the sagittal midline.
            ToSAxisDir = [cosd(ToSangle) 0 sind(ToSangle)] ;
            ToSsagLine = createLine3d(origin,origin+ToSAxisDir) ;
            antToSLine = createLine3d(pBS,teethPt) ;
            [~,antTOSpt,pt2] = distanceLines3d(ToSsagLine,antToSLine) ;

            % calculate the relative position of the resultant force on the
            % sagittal line of the TOS.
            sagTOSLine = createLine3d(origin,antTOSpt) ;
            relPosMuscleOnTOS = (line3dPosition(ptOnToSplane,sagTOSLine)-1).*100  ;
            
            % Caclculate the forces using individual muscle forces
            [Fb,JRFws,JRFbs,M] = biteForce(mResultantForce, ...
                mResultantPos,teethPt,rotTMJ) ;

            % determine the total reaction forces as percentage of the
            % total muscle force
            sumJRFws = sum(JRFws)./vectorNorm3d(totalForce).*100 ;
            sumJRFbs = sum(JRFbs)./vectorNorm3d(totalForce).*100 ;

            % determine the direction of the reaction force
            angle_Fws = atan2d(sumJRFws(3),sumJRFws(1)) ;
            angle_Fbs = atan2d(sumJRFbs(3),sumJRFbs(1)) ;

            % calculate the forces using just the a single vector per
            % muscle
            [FbS,JRFwsS,JRFbsS,Ms] = biteForce(mResultantForce_single, ...
                mResultantPos_single,teethPt,rotTMJ) ;

            % determine the total reaction forces as percentage of the
            % total muscle force
            totalForceSingle = sum(mResultantForce_single);
            sumJRFwsS = sum(JRFwsS)./vectorNorm3d(totalForceSingle).*100 ;
            sumJRFbsS = sum(JRFbsS)./vectorNorm3d(totalForceSingle).*100 ;

            % calculate the difference in reaction forces between using all
            % muscle fibers vs just the middle fiber
            singleMuscleEffect = (sumJRFwsS) - (sumJRFws) ;

            % create the output file (muscleForceDir wsRF bsRF TOSpos effectSingleMuscle )
            outputMuscles(w,1:8,q) = [resForceAngle ...
                sumJRFws([1 3]) sumJRFbs([1 3])...
                relPosMuscleOnTOS singleMuscleEffect([1 3])] ;
        end
    end

    % add the muscle force (PCSA) effect
    effectPCSA = outputMuscles(:,3,q) - outputMuscles(idx_equalMuscles,3,q) ;
    outputMuscles(:,9,q) = effectPCSA ;

    effectMuscle = [outputMuscles(idx_maxMass,3,q) - outputMuscles(idx_equalMuscles,3,q) ...
        outputMuscles(idx_maxMPte,3,q) - outputMuscles(idx_equalMuscles,3,q) ...
        outputMuscles(idx_maxTemp,3,q) - outputMuscles(idx_equalMuscles,3,q)] ;
    muscleEffects(q,:) = effectMuscle ;

    %---------------------------------------------------------------------%
    % PLOT THE RESULTS
    %---------------------------------------------------------------------%

    % PLOT THE MANDIBLE AND POSITION OF WS_RF AND RMF --------------------%

    % For plotting, we are going to plot the forces at different angles
    % with the mandible fixed. 

    invR =(inv(RM)) ; % inverted rotation matrix 

    if q==1
        % plot the mandible
        trisurf(boneModels.Mandible,'edgecolor','none','facecolor',[.8 .8 .8]) ;
        hold on
        axis('vis3d','equal','tight')
        view(34,26)
        rotate3d('on')
        axis off


        % get the limits of the plot
        axisLims = axis ;
        yFactor = (axisLims(4)-axisLims(3))/3 ;

        % plot the trinagle of support
        a3 = trisurf(rotatePatch(TOSpts,invR),'facecolor','r','facealpha',0.3) ;

        % % plot the muscle positions
        % a1 = arrow3(rotateT(mPos,invR), ...
        %     rotateT(mPos+normalizeVector3d(mForce).*yFactor,invR), ...
        %     'f-1',1,2) ;
    end

    % plot the resultant muscle force on the ToS plane. We'll plot the RMF
    % of the equal muscle condition
    set(gca,'ColorOrder',cm(q,:))
    fPos = allMusclePosOnToS(idx_equalMuscles,:) ; 
    fDir = allMuscleDir(idx_equalMuscles,:) ; 
    a2 = arrow3(rotateT(fPos,invR), ...
        rotateT(fPos + fDir.*yFactor,invR),'o-1.5',1,2,1) ;

    % plot the working-side reaction forces
    wsRF = -sumJRFws./100.*yFactor.*3 ;
    p0 = rotateT(pWS,invR) ;

    if wsRF(3)>0
        a5 = arrow3(p0,p0 + wsRF,'o-1.5',1,2) ;
    else
        a5 = arrow3(p0 - wsRF,p0,'o-1.5',1,2) ;
    end
end

camzoom(gca,1.3)
makeHeadlight
axis('auto')

% PLOT HISTOGRAMS OF THE WS RJF AND POSITION OF THE RMF ------------------%
set(figure(102),'windowstyle','docked','color','w') ; clf

p1 = panel() ;
p1.pack(nAngles,4)
p1.de.margin = 5 ;

% determine the limits for the WS reaction force
xLimsJRFz = prctile(outputMuscles(:,3,:),[1 99],"all")' ;

% determine the limits for the muscle position on ToS 
xLimsPos = prctile(outputMuscles(:,6,:),[1 99],"all")' ;

% determine the limits for the position plot
tmp = [xLimsJRFz;xLimsPos] ;
xLimForcePos = [min(tmp(:,1)) max(tmp(:,2))].*1.1 ;

bgColor = ones(1,3).*1 ;
xColor = 'k' ;
yColor = 'w' ;
edgeColor = ones(1,3).*0.1 ;
edgeAlpha = 0.5 ;

nBins = 15 ;
for q = 1:nAngles
    
    % plot the distributions of the vertical ws JRF
    p1(q,2).select() ;
    d1 = outputMuscles(:,3,q) ;
    isBad = d1<xLimsJRFz(1) | d1>xLimsJRFz(2) ;
    d1(isBad) = [] ;
    histogram(d1,nBins,'FaceColor',cm(q,:),'EdgeColor','none' ...
        ,'EdgeAlpha',edgeAlpha) ;
    hold on
    histogram(d1,nBins,'EdgeColor',edgeColor, ...
        'DisplayStyle','stairs') ;
    yLims = ylim ;
    axis([xLimForcePos yLims])
    plot([0 0],yLims,':r','LineWidth',2)
    set(gca,'Color',bgColor,'XColor',xColor,'YColor',yColor,'YTick',[])
    if q~=nAngles
        set(gca,'XTick',[])
    else
        xlabel({'ws JRF','(% total muscle force)'})
    end
    text(xLimForcePos(1)*1.1,yLims(2)*0.5, ...
        [num2str(gapeAngles(q)) ' deg'],"Color",'k', ...
        'HorizontalAlignment','right')
    if q==1
        title('WS Joint Force','color','k')
    end

    % plot the distributions of the resultant muscle force on ToS
    p1(q,3).select() ;
    d2 = outputMuscles(:,6,q) ;
    isBad = d2<xLimsPos(1) | d2>xLimsPos(2) ;
    d2(isBad) = [] ;
    histogram(d2,nBins,'FaceColor',cm(q,:),'EdgeColor','none' ...
        ,'EdgeAlpha',edgeAlpha) ;
    hold on
    histogram(d2,nBins,'EdgeColor',edgeColor, ...
        'DisplayStyle','stairs') ;
    yLims = ylim ;
    axis([xLimForcePos yLims])
    plot([0 0],yLims,':r','LineWidth',2)
    set(gca,'Color',bgColor,'XColor',xColor,'YColor',yColor,'YTick',[])
    if q~=nAngles
        set(gca,'XTick',[])
    else
        xlabel({'RMF position','(% TOS length - 100)'})
    end
    if q==1
        title({'Position RMF','on ToS plane'},'color','k')
    end

end
p1.marginbottom = 18 ;
p1.margintop = 10 ;
p1.fontsize = 12 ;

% PLOT HISROGRAMS OF THE MUSCLE FORCE AND SINGLE FORCE EFFECTS -----------%

% determine the limits for single muscle effect
maxSingle = prctile(abs(outputMuscles(:,8,:)),99,"all") ;

% determine the limits for the effect of muscle force
maxMforce = prctile(abs(outputMuscles(:,9,:)),99,"all") ;

% determine the limits of the effect plots
maxV = max([maxSingle maxMforce]) ;
xLimsEffect = [-maxV maxV].*1.2 ;

% determine the limits for the muscle force angle
xLimsMuscleAngle = [min(outputMuscles(:,1,:),[],'all') ... 
    max(outputMuscles(:,1,:),[],'all')].*1.1 ;

muscleAngleEffect = outputMuscles(:,1,:) - outputMuscles(idx_equalMuscles,1,:) ;
xLimsMuscleAngleEffect = [min(muscleAngleEffect,[],'all') ... 
    max(muscleAngleEffect,[],'all')].*1.1 ;


set(figure(103),'windowstyle','docked','Color','w')
clf
p2 = panel() ;
p2.pack(nAngles,4)
p2.de.margin = 5 ;

for q = 1:nAngles
    % plot the effect of muscle force
    p2(q,2).select() ;
    d1 = outputMuscles(:,9,q) ;
    histogram(d1,nBins,'FaceColor',cm(q,:),'EdgeColor','none' ...
        ,'EdgeAlpha',edgeAlpha) ;
    hold on
    histogram(d1,nBins,'EdgeColor',edgeColor, ...
        'DisplayStyle','stairs') ;
    yLims = ylim ;
    axis([xLimsEffect yLims])
    plot([0 0],yLims,':r','LineWidth',2)
    set(gca,'Color',bgColor,'XColor',xColor,'YColor',yColor,'YTick',[])
    if q~=nAngles
        set(gca,'XTick',[])
    else
        xlabel({'Change in ws JRF','(% total muscle force)'})
    end
    if q==1
        title('Effect Muscle Force','color','k')
    end
    text(xLimsEffect(1)*1.1,yLims(2)*0.5, ...
        [num2str(gapeAngles(q)) ' deg'],"Color",'k', ...
        'HorizontalAlignment','right')

    % plot the effect of single muscle use
    p2(q,3).select() ;
    d1 = outputMuscles(:,8,q) ;
    histogram(d1,nBins,'FaceColor',cm(q,:),'EdgeColor','none' ...
        ,'EdgeAlpha',edgeAlpha) ;
    hold on
    histogram(d1,nBins,'EdgeColor',edgeColor, ...
        'DisplayStyle','stairs') ;
    yLims = ylim ;
    axis([xLimsEffect yLims])
    plot([0 0],yLims,':r','LineWidth',2)
    set(gca,'Color',bgColor,'XColor',xColor,'YColor',yColor,'YTick',[])
    if q~=nAngles
        set(gca,'XTick',[])
    else
        xlabel({'Change in ws JRF','(% total muscle force)'})
    end
    if q==1
        title('Effect Single Muscle','color','k')
    end
    
end
p2.marginbottom = 18 ;
p2.margintop = 10 ;
p2.fontsize = 12 ;

% save the figure and results
prts = strsplit(fileInfo1.name,'_') ;
expFilename = ['TOSmuscle_' strjoin(prts(2:4),'_') '.png'] ;
expFilename2 = ['TOSmuscle2_' strjoin(prts(2:4),'_') '.png'] ;
% exportgraphics(figure(101),expFilename,"Resolution",300,"BackgroundColor","current")
% exportgraphics(figure(102),expFilename2,"Resolution",300,"BackgroundColor","current")
% save('TOS_results.mat','outputMuscles','muscleEffects','-append')