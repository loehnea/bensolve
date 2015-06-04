% Example: MOLP with q=5, unbounded,
% recession cone of upper image has 22 extreme directions (main effort in phase 1)

clear('vlp');

vlp.B=[1 1 1 1 1; 2 1 1 1 1; 1 2 1 1 1; 1 1 2 1 1; 1 1 1 2 1; 1 1 1 1 2; 2 2 1 1 1; 2 1 2 1 1; 2 1 1 2 1; 2 1 1 1 2; 1 2 2 1 1; 1 2 1 2 1; 1 2 1 1 2; 1 1 2 2 1; 1 1 2 1 2; 1 1 1 2 2; 2 2 2 1 1; 2 2 1 2 1; 2 2 1 1 2; 2 1 2 1 2; 2 1 1 2 2; 1 2 2 2 1; 1 2 1 2 2; 1 2 2 1 2; 1 2 2 2 1; 1 1 2 2 2; 1 2 2 2 2; 2 1 2 2 2; 2 2 1 2 2; 2 2 2 1 2; 2 2 2 2 1 ];
vlp.a=[1; zeros(30, 1)];
vlp.P=eye(5,5);

prob2vlp(vlp,'ex11.vlp');
