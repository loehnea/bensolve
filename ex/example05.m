% Example: VLP with q = 3 and 4 generating vectors of C

% see http://bensolve.org/demo.html

clear('vlp');

vlp.B=[ones(1,3);1 2 2;2 2 1;2 1 2];	% coefficient matrix
vlp.a=[1;3/2;3/2;3/2];					% lhs of constraints
vlp.P=[1 0 1; 1 1 0; 0 1 1];			% objective matrix
vlp.l=[0;0;0];							% lower variable bounds

% generating vectors of ordering cone C
 vlp.Y=[1 0 0; 0 1 0; -1 0 2; 0 -1 2]';

% alternative variant: generating vectors of the dual of the ordering cone C
% vlp.Z=[2 2 1 ; 2 0 1 ; 0 0 1 ; 0 2 1]';

% duality parameter vector (must belong to interior of C)
% if not given, it is computed automatically
% dual problem depends on c
vlp.c=[1;1;1];

prob2vlp(vlp,'ex05.vlp');
