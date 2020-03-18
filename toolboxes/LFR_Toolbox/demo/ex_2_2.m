% Linear systems as SS-objects, PCK-objects and constant matrices
% can be converted to LFR-objects by invoking the function "convert2lfr"

% The uncertainty blocks are defined as empty ones except the dynamic one.

   sys0 = rss(5,2,3);
   [a,b,c,d] = ssdata(sys0);

   sys1 = lfr(a,b,c,d,'c');
   sys2 = lfr(sys0);
   sys3 = lfr(pck(a,b,c,d));

% The three above systems are strictly equivalent.

% Conversion between LFR-objects and UMAT-objects.
% UMAT-object:

   x=ureal('a',2.2,'Range',[1 3]);
   y=ureal('b',3.3,'Range',[2 4]);
   z=ureal('c',5.5,'Range',[3 7]);
   s=ucomplex('d',1,'Radius',2);
   t=ucomplex('e',i,'Radius',3);
   A=[x*y+2 z*y;s 1+t*z]

% Conversion to LFR-object using the function lfr

   B = lfr(A);
   size(B)

% Back conversion to UMAT-object

   C = lfr2rob(B);
   distlfr(A,C)
