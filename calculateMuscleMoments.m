%% calculateMuscleMoments

% This script calculates the moment arms and relative lengths of the muscle
% segments of the three jaw addutor muscles at different gape angles.
%
% The script requires:
%   BONE MODEL and ROTATION DATA:
%       mat files that contains the 3D triangulated cranium and mandible
%       models and the rotation matrices for each gape angle. These data
%       are exported from the 'modelPreparation' MATLAB app
%   MUSCLE DATA:
%       mat file that contains the position of the three jaw adductor
%       muscles at occlusion (0-degree gape). This file is obtained from
%       the 'muscleSelection' MATLAB app
% 
% The script exports a 'MuscleMoments_*.mat' file that contains the position,
% relative lengths, and moment arms of the jaw adductor muscles.

clear
clc

% let's start by loading the MUSCLE DATA
% The required file comes from 'muscleSelection' app
muscleFile = dir('muscleData_*.mat') ;
load(muscleFile.name)

% Rotation and muscle model data come from the 'stlModelPreparation' app

% now let's load the ROTATION DATA
rotationFile = dir('rotations_*.mat') ;
% rotationFile = dir('RotationsData_*.mat') ;
load(rotationFile.name) ;

% now let's load the BONE MODELS
boneFile = dir('boneModels_*.mat') ;
load(boneFile.name) ;

%% get and prepare the bone and muscle data
% prepare the muscle data
muscles = muscleData.muscles ;
muscleNames = fieldnames(muscles) ;

modelCranium  = cleanMesh(boneModels.Cranium) ;
modelMandible = cleanMesh(boneModels.Mandible) ;

%% calculate the muscle position and direction for each rotation
nRotations = size(rotations,1) ;

muscleMomentData = [] ;
nPoints = 20 ; % number of points per muscle segment

muscleMoment = table(rotations(:,1),rotations(:,2),rotations(:,3), ...
    'VariableNames',{'GapeAngle','xTransl','zTransl'}) ;
musclesPerAngle = cell(11,3) ;
for i = 1:nRotations
    multiWaitbar('Gape Angle',i/nRotations)

    % get the rotation for that gape
    rotValues = rotations(i,:) ;
    Ry = troty(rotValues(1),'deg') ;
    T  = transl(rotValues(2),0,rotValues(3)) ;
    T_Rot = T*Ry ;

    % rotate the mandible model and create a combined variable with both
    % models
    rotModelMandible = rotatePatch(modelMandible,T_Rot) ;
    combinedModels.Cranium  = modelCranium ;
    combinedModels.Mandible = rotModelMandible ;

    % get the position of the temporal fossa
    if i==1
        % see if the temporalis is on the left or right side of the model
        tempAttchment  = reshape(muscles.Temporalis(end,:),3,[])' ;
        if tempAttchment(1,2)<0; isRight = true ; else; isRight = false ; end

        % get the the points that define the edge of the temporal fossa
        fossaPts = getTemporalFossaEdge(rotModelMandible,modelCranium, ...
            isRight) ;
    end

    allMuscles = cell(1,3) ;
    for j = 1:3 % for the different muscles
        muscleName = muscleNames{j} ;
        musclePos  = muscles.(muscleName) ;
        nSegments = size(musclePos,2)/3 ;

        if strcmp(muscleName,'Temporalis')
            isTemporalis = true ;
        else
            isTemporalis = false ;
        end
    
        % get the cranium and mandible attachments
        craniumPts  = reshape(musclePos(1,:),3,[])' ;
        mandiblePts = rotateT(reshape(musclePos(end,:),3,[])',T_Rot) ; 

        % initialize rotated muscle variables
        musclePts = nan(size(musclePos)) ;
        
        for q = 1:nSegments
            multiWaitbar('Muscle segments',((j-1)*nSegments + nSegments)/(3*nSegments))

            % calculate the muscle segment
            muscleSegmentPts = getMuscleSegments(combinedModels,...
                craniumPts(q,:),mandiblePts(q,:),nPoints,muscleName,fossaPts) ;
            
            musclePts(:,3*q-2:3*q) = muscleSegmentPts ;
        end

        allMuscles{j} = musclePts ;
    end

    % make sure that the points are on the surface of the model
    % Combine the models
    FV_md = triang2mesh(combinedModels.Mandible) ;
    FV_cr = triang2mesh(combinedModels.Cranium) ;
    [v,f] = concatenateMeshes(FV_cr,FV_md) ;
    TR = triangulation(f,v) ;

    % combine the muscles
    allMusclePts = cell2mat(allMuscles);
    sizeRef = size(allMusclePts) ;
    allMusclePts = reshape(allMusclePts',3,[])' ;
    [distPts,projPts] = closestPoint2Mesh(TR,allMusclePts) ;

    % check the average size of the faces
    meanFaceSize = faceSize(TR.ConnectivityList,TR.Points) ;
    critDist = meanFaceSize/10 ;

    isOnSurface = distPts<critDist ;
    finalPts = nan(size(allMusclePts)) ;
    finalPts( isOnSurface,:) = projPts(   isOnSurface,:) ;
    finalPts(~isOnSurface,:) = allMusclePts(~isOnSurface,:) ;
    
    outputMuscle = reshape(finalPts',[],nPoints)' ;

    nF = nSegments*3 ;
    finalMuscles = mat2cell(outputMuscle,nPoints,[nF nF nF]) ;
    musclesPerAngle(i,:) = finalMuscles ;
end
multiWaitbar('closeall')
muscleMoment.musclePosition = musclesPerAngle ;
muscleMoment.Properties.VariableDescriptions = ["","","",...
    "Masseter - MedPterygoid - Temporalis"] ;

%% calculate the muscle moments

% initialize the output variables
Moments = cell(nRotations,3) ;
Lengths = cell(nRotations,3) ;

for i = 1:nRotations
    for j = 1:3
        % get the muscle position data
        musclePoints = muscleMoment.musclePosition{i,j} ;

        % calculate the muscle length and the moments for each rotation
        Origin = zeros(1,3) ;
        [mMoments,  mLength] = muscleLenghtMoment(musclePoints,Origin) ;
        Moments{i,j} = mMoments ;
        Lengths{i,j} = mLength ;
    end
end
muscleMoment.muscleLength = Lengths ;
muscleMoment.muscleMoment = Moments ;

% save the data
path1 = muscleFile.folder ;
[~,species] = fileparts(path1) ;
save(['MuscleMoments_' species '.mat'],'muscleMoment')

%% visualize the results
resolution = get(0,"ScreenSize") ;
maxX = resolution(3)*0.95 ; maxY = resolution(4)*0.8 ;
posX = resolution(3)*0.02 ; posY = resolution(4)*0.05 ;
set(figure,'color','k','Position',[posX,posY,maxX,maxY]) ;
nRotations = size(rotations,1) ;

% check if the muscles are on the left or the right side. If they are on
% the right side, invertD the Y axis so we can see them.
massLoc = sign(mean(muscleMoment.musclePosition{1,1}(:,2))) ;
pterLoc = sign(mean(muscleMoment.musclePosition{1,2}(:,2))) ;
tempLoc = sign(mean(muscleMoment.musclePosition{1,3}(:,2))) ;
if mean([massLoc pterLoc tempLoc]) > 0 % is on the right side
    doInverseY = true ;
else
    doInverseY = false ;
end

idx = [1:5 7:11] ;
for w = 2:nRotations
    subplot(2,6,idx(w-1))

    rotValues = rotations(w,:) ;
    Ry = troty(rotValues(1),'deg') ;
    T  = transl(rotValues(2),0,rotValues(3)) ;
    T_Rot = T*Ry ;

    trisurf(modelCranium,'facecolor',[.8 .8 .8],'edgecolor','none') ; hold on
    trisurf(rotatePatch(modelMandible,T_Rot),'facecolor',[.8 .8 .8],...
        'edgecolor','none') ;
    view(0,0)
    axis("vis3d",'equal','tight', 'off')

    p3(muscleMoment.musclePosition{w,1},'parula','o-k','markersize',5,'linewidth',1.5) ;
    p3(muscleMoment.musclePosition{w,2},'parula','o-k','markersize',5,'linewidth',1.5) ;
    p3(muscleMoment.musclePosition{w,3},'parula','o-k','markersize',5,'linewidth',1.5) ;

    if doInverseY
        set(gca,'Ydir','reverse')
    end
end
makeHeadlight

% plot the muscle length vs gape
cm = colormap(parula(7)) ;
subplot(3,6,6) ;
colororder(cm)
L = cell2mat(muscleMoment.muscleLength(:,1)) ;
plot(rotations(:,1),L./L(1,:).*100,'LineWidth',2)
title('Masseter','Color','w')
ylabel('Muscle Length (% L_0)')
ax1 = gca ;
ax1.Color = 'k' ;
ax1.YColor = 'w' ;
ax1.XColor = 'w' ;
ax1.FontSize = 12 ;

subplot(3,6,12)
colororder(cm)
L = cell2mat(muscleMoment.muscleLength(:,2)) ;
plot(rotations(:,1),L./L(1,:).*100,'LineWidth',2)
title('Med Pterygoid','Color','w')
ylabel('Muscle Length (% L_0)')
ax2 = gca ;
ax2.Color = 'k' ;
ax2.YColor = 'w' ;
ax2.XColor = 'w' ;
ax2.FontSize = 12 ;

subplot(3,6,18)
colororder(cm)
L = cell2mat(muscleMoment.muscleLength(:,3)) ;
plot(rotations(:,1),L./L(1,:).*100,'LineWidth',2)
title('Temporalis','Color','w')
ylabel('Muscle Length (% L_0)')
xlabel('Gape (degrees)')
ax3 = gca ;
ax3.Color = 'k' ;
ax3.YColor = 'w' ;
ax3.XColor = 'w' ;
ax3.FontSize = 12 ;

% exportgraphics(gcf,'Rotations muscle.png','BackgroundColor','current')