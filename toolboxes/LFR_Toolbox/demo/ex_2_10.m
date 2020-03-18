% Well-posedness and non-singularity radii

   lfrs a b c d
   M = [(1+a)/(2-b-c) 2;2*a 3+d];

% M is not well-posed for b + c = 2, worst case for b = c = 1

   [rad_min,rad_max,pertnames,pert] = wp_rad(M)

% M becomes singular for a = 1/3 and b = c = d = -1/3

   [rad_min,rad_max,pertnames,pert] = ns_rad(M)
