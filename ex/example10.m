% The 'bensolvehedron', see the titlepage of the reference manual
%
% MOLP with q objectives, n=(q+2*m)^q variables and constraints
%
% To compute the bensolvehedron, type
%   ./bensolve ex/ex10.vlp -p
% Plot the resulting OFF files (or INST files), for instance, with GEOMVIEW:
%    geomview ex/ex10_p.inst   OR    geomview ex/ex10_p.off (unscaled version)
%    geomview ex/ex10_d.inst   OR    geomview ex/ex10_d.off (unscaled version)
% add green color (!)

% adopt q and m to generate other (larger) problems
q=3;
m=2;

n=(q+2*m)^q;

clear('vlp');

% feasible set is n dimensional hyper cube
vlp.B=eye(n);
vlp.a=zeros(n,1);
vlp.b=ones(n,1);

% objective map
P=zeros(n,q);
for i=1:n
	line = dec2base(i-1,q+2*m,q)-'0';
	line = line - (q+2*m-1)/2;
	P(i,:) = line;
end
vlp.P = P';

prob2vlp(vlp,'ex10.vlp');