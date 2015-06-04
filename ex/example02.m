% Example: MOLP with 2 objectives which is infeasible

% v-min [x1;x2]
%
% 0 <= 3*x1 +   x2 <= 1
% 0 <=   x1 + 2*x2 <= 1
% 1 <=   x1 +   x2 <= 2

clear('vlp');

vlp.B=[[3 1];[1 2];[1 1]];
vlp.b=[1; 1; 2];
vlp.a=[0; 0; 1];
vlp.P=[1 0; 0 1];

prob2vlp(vlp,'ex02.vlp');
