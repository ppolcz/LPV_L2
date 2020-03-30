function grd=multigrid(varargin)
%It does the same job as ndgrid, but returns the grid
%points in the rows of grd. The number of inputs is free.
%Syntax: grd=MULTIGRID(...)
%Example:
%   grd = multigrid([-1,2],[3 4],[5 6])
%   grd =       -1     3     5
%                2     3     5
%               -1     4     5
%                2     4     5
%               -1     3     6
%                2     3     6
%               -1     4     6
%                2     4     6

d=length(varargin);
varargin=flip(varargin);
blks=cell(1,d);
[blks{:}]=ndgrid(varargin{:});
ne=numel(blks{1});
grd=zeros(ne,d);
for k=1:d
    grd(:,k)=reshape(blks{k},ne,1);
end;
grd=flip(grd,2);


