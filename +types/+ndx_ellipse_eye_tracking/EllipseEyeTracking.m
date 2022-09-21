classdef EllipseEyeTracking < types.core.EyeTracking & types.untyped.GroupClass
% ELLIPSEEYETRACKING Stores detailed eye tracking information output from DeepLabCut


% PROPERTIES
properties
    corneal_reflection_tracking; % corneal reflection tracking
    eye_tracking; % eye tracking
    likely_blink; % Indicator of whether there was a probable blink for this frame
    pupil_tracking; % pupil tracking
end

methods
    function obj = EllipseEyeTracking(varargin)
        % ELLIPSEEYETRACKING Constructor for EllipseEyeTracking
        %     obj = ELLIPSEEYETRACKING(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % corneal_reflection_tracking = EllipseSeries
        % eye_tracking = EllipseSeries
        % likely_blink = TimeSeries
        % pupil_tracking = EllipseSeries
        obj = obj@types.core.EyeTracking(varargin{:});
        
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'corneal_reflection_tracking',[]);
        addParameter(p, 'eye_tracking',[]);
        addParameter(p, 'likely_blink',[]);
        addParameter(p, 'pupil_tracking',[]);
        misc.parseSkipInvalidName(p, varargin);
        obj.corneal_reflection_tracking = p.Results.corneal_reflection_tracking;
        obj.eye_tracking = p.Results.eye_tracking;
        obj.likely_blink = p.Results.likely_blink;
        obj.pupil_tracking = p.Results.pupil_tracking;
        if strcmp(class(obj), 'types.ndx_ellipse_eye_tracking.EllipseEyeTracking')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.corneal_reflection_tracking(obj, val)
        obj.corneal_reflection_tracking = obj.validate_corneal_reflection_tracking(val);
    end
    function obj = set.eye_tracking(obj, val)
        obj.eye_tracking = obj.validate_eye_tracking(val);
    end
    function obj = set.likely_blink(obj, val)
        obj.likely_blink = obj.validate_likely_blink(val);
    end
    function obj = set.pupil_tracking(obj, val)
        obj.pupil_tracking = obj.validate_pupil_tracking(val);
    end
    %% VALIDATORS
    
    function val = validate_corneal_reflection_tracking(obj, val)
        val = types.util.checkDtype('corneal_reflection_tracking', 'types.ndx_ellipse_eye_tracking.EllipseSeries', val);
    end
    function val = validate_eye_tracking(obj, val)
        val = types.util.checkDtype('eye_tracking', 'types.ndx_ellipse_eye_tracking.EllipseSeries', val);
    end
    function val = validate_likely_blink(obj, val)
        val = types.util.checkDtype('likely_blink', 'types.core.TimeSeries', val);
    end
    function val = validate_pupil_tracking(obj, val)
        val = types.util.checkDtype('pupil_tracking', 'types.ndx_ellipse_eye_tracking.EllipseSeries', val);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.EyeTracking(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        if ~isempty(obj.corneal_reflection_tracking)
            refs = obj.corneal_reflection_tracking.export(fid, [fullpath '/corneal_reflection_tracking'], refs);
        else
            error('Property `corneal_reflection_tracking` is required in `%s`.', fullpath);
        end
        if ~isempty(obj.eye_tracking)
            refs = obj.eye_tracking.export(fid, [fullpath '/eye_tracking'], refs);
        else
            error('Property `eye_tracking` is required in `%s`.', fullpath);
        end
        if ~isempty(obj.likely_blink)
            refs = obj.likely_blink.export(fid, [fullpath '/likely_blink'], refs);
        else
            error('Property `likely_blink` is required in `%s`.', fullpath);
        end
        if ~isempty(obj.pupil_tracking)
            refs = obj.pupil_tracking.export(fid, [fullpath '/pupil_tracking'], refs);
        else
            error('Property `pupil_tracking` is required in `%s`.', fullpath);
        end
    end
end

end