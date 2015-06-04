% Example: MOLP with 2 objectives the upper image of which has no vertex

% v-min [x1;x2]
%
% 1 <= x1 + x2 + x3
% 1 <= x1 + x2 - x3

clear('vlp');

vlp.B=[[1 1 1];[1 1 -1]];
vlp.a=[1;1];
vlp.P=[1 0 0; 0 1 0];

prob2vlp(vlp,'ex03.vlp');
