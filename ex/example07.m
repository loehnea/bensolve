% Example: MOLP with 3 objectives, 1211 constraints and 1143 variables
%
% Example (PL) in
% 
% Shao, L., Ehrgott, M.: Approximately solving multiobjective linear programmes in objective space and
% an application in radiotherapy treatment planning. Math. Methods Oper. Res. 68(2), 257Ð276 (2008)
%
% enlarge epsilon in phase 2 and use primal simplex algorithm, for instance, run
% ./bensolve ex/ex07.vlp -m 2 -e 0.05 -l primal_simplex

clear('vlp');

vlp.B=load('example07.txt');
vlp.P = -[zeros(1,1140),-1,0,0;zeros(1,1140),0,-1,0;zeros(1,1140),0,0,-1];
vlp.l = [zeros(1140,1);0;-45;0];
vlp.s = [Inf*ones(1140,1);17.07;12;90.64]; 
vlp.b =[90.64*ones(67,1);-85.3601*ones(67,1);60*ones(37,1);45*ones(8,1);60*ones(46,1);ones(986,1)];

prob2vlp(vlp,'ex07.vlp');
