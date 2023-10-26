function h = cube(limits, varargin)
% CUBE Plot a cube
%   h = CUBE(limits) plots a cube. limits is a 2 x 3 matrix specifying the
%       limits (min, max) in x, y and z for the cube.
%   h = CUBE(limits, ...) additional arguments accepted by plot3 can be
%       passed, such as line style, marker symbol and color.
    
    x = repelem(limits(:,1),4,1);
    y = repelem(limits(:,2),2,1);
    y = [y; y];
    z = repmat(limits(:,3),4,1);
    
    inds = [...
        1 1 1 2 2 3 3 4 5 5 6 7;...
        2 3 5 4 6 4 7 8 6 7 8 8];
    
    h = plot3(x(inds), y(inds), z(inds), varargin{:});
end

