% PLOTLFR      - 1-D or 2-D plot of an LFR-object entry
%-----------------------------------------------------------------
% PURPOSE
% For  1-D  or  2-D  plotting  of one entry (real case) or modulus
% (complex case) of a non-dynamic lfr-objects.
%
% SYNOPSIS
% plotlfr(lfrobj,px[,py][,lfrind])
%
% DESCRIPTION
% SISO  case:  plotlfr plots the values of the lfr-object lfrobj.
% MIMO  case:  plotlfr  plots  the values of one entry of lfrobj,
% in this case, this entry is defined by lfrind.
% In  the 1-D case, px defines the x-axis (parameter to be gridded
% and variation bounds).  In the 2-D case,  py defines the y-axis.
% If  the  plotted  object  depends  on  more  than  2 parameters,
% additional  parameter  values are fixed to zero (can be fixed to
% other values before invoking plolfr).
%
% INPUT ARGUMENTS
% lfrobj lfr-object.
% px     Cell of the form {name1,min1,max1,npt1} in which
%        * name1 (string  or  lfr-object) indicates the uncertain
%          parameter that varies on the x-axis,
%        * [min1 max1] is the interval of variation,
%        * npt1 is the  number  of points in this interval.
% py     Similar to px.
% lfrind 1_by_2  integer  matrix the entries of which define which
%        entry of lfrobj is plotted (MIMO case).
%#-----------------------------------------------------------------
% % EXAMPLE
%    lfrs a b c
%    lfrobj = [1/(1+a) b^2 a+b+c;1+a^2 c*2*j/(1+a+b) 4];
%
% % lfrobj(1,1) = 1/(1+a)
%    plotlfr(lfrobj,{a,-0.9,4,20},[1 1])
%    plotlfr(lfrobj,{'a',-0.9,4,20},[1 1])
%
% % lfrobj(1,3) = a + b +c at c = 3;
%    c = 3; lfrobj2 = eval(lfrobj);
%    plotlfr(lfrobj2,{a,0,4,10},{b,0,4,10},[1 3])
%
% % Modulus of lfrobj(2,2) at b = 0 (default value)
%    plotlfr(lfrobj,{'a',0,4,10},{'c',0,4,10},[2 2])
