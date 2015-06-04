% VLP with 3 objectives, 4608 constrains and 36939 variables
%
% Example 6.6 in
% Hamel, A.H.; Loehne,A.; Rudloff,B.: A Benson-type algorithm for linear vector
% optimization and applications. J. Glob. Optim. 59, No. 4, 811-836 (2014)
%
% special options are necessary to run example:
%
% For instance, run
% ./bensolve ex/ex09.vlp -e 1e-2 -m 3 -L primal_simplex -l primal_simplex -p
%
% (-e: set epsilon in Phase 2)
% (-m set message level to have more output on screen)
% (-L use primal simplex algorithm in Phase 1 of Benson's algorithm)
% (-l use primal simplex algorithm in Phase 2 of Benson's algorithm)
% (-p generate graphics files)
%
% alternatively, run
%
% ./bensolve ex/ex09.vlp -e 1e-2 -m3 -A dual -a dual -p
%
% (-A use dual Benson algorithm in Phase 1)
% (-a use dual Benson algorithm in Phase 2)
% (-p generate graphics files)

clear('vlp');
load('example09');
vlp.c=[1;1;1];
prob2vlp(vlp,'ex09.vlp');