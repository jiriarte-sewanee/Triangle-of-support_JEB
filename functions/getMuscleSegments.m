function musclePts = getMuscleSegments(model,...
    craniumPts,mandiblePts,nPtsPerSegment,muscle,fossaPts)

% GETMUSCLESEGMENTS calculates the muscles segments that connect the
% CRANIUMPTS and the MANDIBLEPTS. If the axis connecting the cranium and
% mandible points passes through the bone, the muscle segment wraps around
% the model. NPTSPERSEGMENT indicates how many points would form each
% segment. MUSCLE is a character array that indicates which muscle
% (masseter, medPterygoid, or temporalis) is being preocessed. FOSSAPTS is
% a series of points that indicate the posterior edge of the temporal
% fossa. This is needed to make that the temporalis muscle passes through
% the fossa, especially at large gaopes.

if size(craniumPts)~=size(mandiblePts)
    error('The mandible and cranium should have the same number of points')
end

if size(craniumPts,1)==1; craniumPts = reshape(craniumPts',3,[])' ; end
if size(mandiblePts,1)==1; mandiblePts = reshape(mandiblePts',3,[])' ; end
    
% get the number of muscle segments to analyze
nMuscleSegments = size(craniumPts,1) ;

% initialize output variable
musclePts = nan(nPtsPerSegment,nMuscleSegments*3) ;

% define the models
TR_cr = model.Cranium ;
TR_md = model.Mandible ;

for i = 1:nMuscleSegments

    % get the muscle points
    craniumPt = craniumPts (i,:) ;
    mandiblePt = mandiblePts(i,:) ;

    % find the non-convex points on both the cranium and the mandible
    % that intersect the plane defined by the proximal and distal
    % muscle points

    % for the CRANIUM ----------------------------------------------------%
    % define the muscle axis, the line between tha cranium and mandible
    % attachment point
    muscleAxisPts = [craniumPt; mandiblePt] ;
    isCranium = true ;
    surfPts_cranium = getMuscleOutline(TR_cr, muscleAxisPts, ...
        muscle,isCranium,fossaPts) ;

    % for the MANDIBLE ---------------------------------------------------%
    % define the muscle axis, the line between tha mandible point and the
    % last calculated cranial point
    muscleAxisPts_md = [mandiblePt; surfPts_cranium(end,:)] ;
    isCranium = false ;
    surfPts_mandible = getMuscleOutline(TR_md, muscleAxisPts_md, ...
        muscle,isCranium,[]) ;

    % Procese the final data ---------------------------------------------%
    % export the calculated surface points
    combinedPts = [surfPts_cranium; flipud(surfPts_mandible)] ;

    interpPts = curvspace(combinedPts,nPtsPerSegment) ;

    musclePts(:,3*i-2:3*i) = interpPts ;
end