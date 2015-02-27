function [t_interior,f_interior] = extract_interior(t_data,f_data,t_min,t_max)
%EXTRACT_INTERIOR - EXTRACTS all points in range tmin<= t_data <= tmax, and
%corresponding points of f (of same length). extract_interior_V4.m is
%identical.
%
% Syntax: [t_interior,f_interior] = extract_interior(t_data,f_data,t_min,t_max)
%
% Inputs:
%    t_data - time array matched to f_data, monotonically increasing
%    f_data - data array matched to t_data
%    t_min - starting value in t_data to extract points from.
%            Same physical units as t_data.
%    t_max - ending value in t_data to extract points from
%
% Outputs:
%    t_interior - Elements of t_data between t_min and t_max
%    f_interior - Elements of f_data matching t_interior
%
% Example: 
%    [t_interior,f_interior] = extract_interior([0 1 2 3],[0 3 6 9],1,3)
%    t_interior equals [1 2 3], f_interior equals [3 6 9].
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Joseph P. Feser
% Work address --
% email: --
% Website: --
% [Month] [Year]; Last revision: unknown

%------------- BEGIN CODE --------------
    Ndata=length(t_data);
    Nmin=Ndata;
    for i=Ndata:-1:1
        if t_data(i)>=t_min
            Nmin=i;
        end
    end
    Nmax=Nmin;
    for i=Nmin:Ndata
        if t_data(i)<=t_max
            Nmax=i;
        end
    end
    t_interior=t_data(Nmin:Nmax,1);
    f_interior=f_data(Nmin:Nmax,1);
end

