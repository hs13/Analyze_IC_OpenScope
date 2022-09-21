classdef EllipseSeries < types.core.SpatialSeries & types.untyped.GroupClass
% ELLIPSESERIES Information about an ellipse moving over time


% PROPERTIES
properties
    angle; % angle that ellipse is rotated by (phi)
    area; % ellipse area, with nan values in likely blink times
    area_raw; % ellipse area, with no regard to likely blink times
    height; % height of ellipse
    width; % width of ellipse
end

methods
    function obj = EllipseSeries(varargin)
        % ELLIPSESERIES Constructor for EllipseSeries
        %     obj = ELLIPSESERIES(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % angle = float
        % area = float
        % area_raw = float
        % height = float
        % width = float
        varargin = [{'data_conversion' types.util.correctType(1.0, 'float32') 'data_resolution' types.util.correctType(-1.0, 'float32') 'data_unit' 'meters'} varargin];
        obj = obj@types.core.SpatialSeries(varargin{:});
        
        
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'angle',[]);
        addParameter(p, 'area',[]);
        addParameter(p, 'area_raw',[]);
        addParameter(p, 'height',[]);
        addParameter(p, 'width',[]);
        misc.parseSkipInvalidName(p, varargin);
        obj.angle = p.Results.angle;
        obj.area = p.Results.area;
        obj.area_raw = p.Results.area_raw;
        obj.height = p.Results.height;
        obj.width = p.Results.width;
        if strcmp(class(obj), 'types.ndx_ellipse_eye_tracking.EllipseSeries')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.angle(obj, val)
        obj.angle = obj.validate_angle(val);
    end
    function obj = set.area(obj, val)
        obj.area = obj.validate_area(val);
    end
    function obj = set.area_raw(obj, val)
        obj.area_raw = obj.validate_area_raw(val);
    end
    function obj = set.height(obj, val)
        obj.height = obj.validate_height(val);
    end
    function obj = set.width(obj, val)
        obj.width = obj.validate_width(val);
    end
    %% VALIDATORS
    
    function val = validate_angle(obj, val)
        val = types.util.checkDtype('angle', 'float', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = val.dims;
        else
            valsz = size(val);
        end
        validshapes = {[1]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_area(obj, val)
        val = types.util.checkDtype('area', 'float', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = val.dims;
        else
            valsz = size(val);
        end
        validshapes = {[1]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_area_raw(obj, val)
        val = types.util.checkDtype('area_raw', 'float', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = val.dims;
        else
            valsz = size(val);
        end
        validshapes = {[1]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_data(obj, val)
        val = types.util.checkDtype('data', 'numeric', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = val.dims;
        else
            valsz = size(val);
        end
        validshapes = {[2 Inf]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_data_continuity(obj, val)
        val = types.util.checkDtype('data_continuity', 'char', val);
    end
    function val = validate_data_conversion(obj, val)
        val = types.util.checkDtype('data_conversion', 'float32', val);
    end
    function val = validate_data_resolution(obj, val)
        val = types.util.checkDtype('data_resolution', 'float32', val);
    end
    function val = validate_data_unit(obj, val)
        val = types.util.checkDtype('data_unit', 'char', val);
    end
    function val = validate_height(obj, val)
        val = types.util.checkDtype('height', 'float', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = val.dims;
        else
            valsz = size(val);
        end
        validshapes = {[1]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_width(obj, val)
        val = types.util.checkDtype('width', 'float', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = val.dims;
        else
            valsz = size(val);
        end
        validshapes = {[1]};
        types.util.checkDims(valsz, validshapes);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.SpatialSeries(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        if ~isempty(obj.angle)
            if startsWith(class(obj.angle), 'types.untyped.')
                refs = obj.angle.export(fid, [fullpath '/angle'], refs);
            elseif ~isempty(obj.angle)
                io.writeDataset(fid, [fullpath '/angle'], obj.angle);
            end
        else
            error('Property `angle` is required in `%s`.', fullpath);
        end
        if ~isempty(obj.area)
            if startsWith(class(obj.area), 'types.untyped.')
                refs = obj.area.export(fid, [fullpath '/area'], refs);
            elseif ~isempty(obj.area)
                io.writeDataset(fid, [fullpath '/area'], obj.area);
            end
        else
            error('Property `area` is required in `%s`.', fullpath);
        end
        if ~isempty(obj.area_raw)
            if startsWith(class(obj.area_raw), 'types.untyped.')
                refs = obj.area_raw.export(fid, [fullpath '/area_raw'], refs);
            elseif ~isempty(obj.area_raw)
                io.writeDataset(fid, [fullpath '/area_raw'], obj.area_raw);
            end
        else
            error('Property `area_raw` is required in `%s`.', fullpath);
        end
        if ~isempty(obj.height)
            if startsWith(class(obj.height), 'types.untyped.')
                refs = obj.height.export(fid, [fullpath '/height'], refs);
            elseif ~isempty(obj.height)
                io.writeDataset(fid, [fullpath '/height'], obj.height);
            end
        else
            error('Property `height` is required in `%s`.', fullpath);
        end
        if ~isempty(obj.width)
            if startsWith(class(obj.width), 'types.untyped.')
                refs = obj.width.export(fid, [fullpath '/width'], refs);
            elseif ~isempty(obj.width)
                io.writeDataset(fid, [fullpath '/width'], obj.width);
            end
        else
            error('Property `width` is required in `%s`.', fullpath);
        end
    end
end

end