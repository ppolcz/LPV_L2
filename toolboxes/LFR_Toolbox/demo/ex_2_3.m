% Use of the function lfr for defining Delta-blocks

% A 4-by-4 scalar repeated block
   X = lfr('x','ltisr',4,[2 8],'minmax');
   size(X)

% A complex scalar uncertainty in the disc of center 1+i and radius 2.2
   Y = lfr('y','ltisc',1,[1+i 2.2],'disc');
   size(Y)

% A full complex block with frequency weighted bounds
   Z = lfr('z','ltifc',[2 4],ltisys(-10,-10,1,0),'freq');
   size(Z)

