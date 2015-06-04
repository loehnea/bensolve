% Example: VLP with 2 objectives

% max [x1 - x2; x1 + x2]
%
% w.r.t. cone C ={(y1,y2) | 2 y1 -y2 >= 0, -y1 + 2 y2 >= 0}
% 
% 1 <= x1 + x2 <= 2 
%
% 0 <= x1 <= 1
% 0 <= x2

clear('vlp');

vlp.opt_dir=-1; 	% maximization
vlp.Z=[2 -1; -1 2];	% generators of dual of ordering cone
vlp.c=[1;1];		% geometric duality parameter vector (belongs to interior of C)

vlp.B=[1 1];		% coefficient matrix
vlp.a=1;			% lhs
vlp.b=2;			% rhs
vlp.l=[0;0];		% lower bounds
vlp.s=[1,Inf];		% upper bounds
vlp.P=[1 -1;1 1];	% objective matrix

prob2vlp(vlp,'ex06.vlp');

