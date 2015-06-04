% Example: MOLP with 2 objectives, simplest example

% min [x1 - x2; x1 + x2]
% 
% 6 <= 2*x1 +   x2
% 6 <=   x1 + 2*x2
%
% x1 >= 0
% x2 >= 0

clear('vlp');

vlp.B=[2 1;1 2];    % coefficient matrix
vlp.a=[6;6];        % lower bounds
vlp.P=[1 -1; 1 1];  % objective matrix
vlp.l=[0;0];        % lower variable bounds

prob2vlp(vlp,'ex01.vlp');
