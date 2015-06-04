% This file is part of BENSOLVE - VLP solver
%
% Copyright (C) 2014-2015 Andreas Löhne and Benjamin Weißing
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program (see the reference manual). If not,
% see <http://www.gnu.org/licenses/>

function [ ] = prob2vlp(vlp,filename)
	if isfield(vlp,'ub')&& ~isfield(vlp,'s')
		vlp.s=vlp.ub;
	end
	if isfield(vlp,'lb')&& ~isfield(vlp,'l')
		vlp.l=vlp.lb;
	end
	if ~isfield(vlp,'opt_dir')
		vlp.opt_dir=1;
	end
	if ~isfield(vlp,'Y')
		vlp.Y=[];
	end
	if ~isfield(vlp,'Z')
		vlp.Z=[];
	end
	if ~isfield(vlp,'c')
		vlp.c=[];
	end
	if ~isfield(vlp,'l')
		vlp.l=[];
	end
	if ~isfield(vlp,'s')
		vlp.s=[];
	end
	if ~isfield(vlp,'a')
		vlp.a=[];
	end
	if ~isfield(vlp,'b')
		vlp.b=[];
	end
	if ~isfield(vlp,'b')
		error('Coeffient matrix vlp.B must be given.');
	end
	if ~isfield(vlp,'P')
		error('Objective matrix vlp.P must be given.');
	end
	[fid,message] = fopen(filename,'w');
	if fid==-1
		error('Cannot open file %s \n%s',filename,message);
	end 
	[m,n]=size(vlp.B);
	[q,p]=size(vlp.P);
	if (n<1)
		error('vlp.B must have at least one column.');
	end
	if (q<1)
		error('vlp.P must have at least one row.');
	end
	if (n~=p)
		error('vlp.B and vlp.P must have the same number of columns.');
	end


	% constraint coefficients
	[A_rows,A_cols,A_vals]=find(sparse(vlp.B));
	k=length(A_rows(:));

	% objective coefficients
	[P_rows,P_cols,P_vals]=find(sparse(vlp.P));
	k1=length(P_rows(:));

	% program line and cone coefficients
	if size(vlp.Y,2)>0
		[K_rows,K_cols,K_vals]=find(sparse([vlp.Y]));
		k2=length(K_rows(:));
		str=sprintf(' cone %d %d',size(vlp.Y,2),k2);
	elseif size(vlp.Z,2)>0
		[K_rows,K_cols,K_vals]=find(sparse([vlp.Z]));
		k2=length(K_rows(:));
		str=sprintf(' dualcone %d %d',size(vlp.Z,2),k2);
	else
		str='';
		k2=0;
	end
	
	
	if vlp.opt_dir==1
		opt_dir_str=sprintf('min');
	elseif vlp.opt_dir==-1
		opt_dir_str=sprintf('max');
	else
		error('invalid entry of vlp.opt_dir: use 1 for minimization and -1 for maximization');
	end
	
	% write 'p', 'a', 'k' to file
	fprintf(fid,'p vlp %s %d %d %d %d %d%s\n',opt_dir_str,m,n,k,q,k1,str);
	for i=1:k
		fprintf(fid,'a %d %d %g\n',A_rows(i),A_cols(i),A_vals(i));
	end
	for i=1:k1
		fprintf(fid,'o %d %d %g\n',P_rows(i),P_cols(i),P_vals(i));
	end	
	for i=1:k2
		fprintf(fid,'k %d %d %g\n',K_rows(i),K_cols(i),K_vals(i));
	end
	
	% duality parameter vector
	if ~isempty(vlp.c)
		if size(vlp.c,1)~=q || size(vlp.c,2) ~= 1
			error('vlp.c has wrong dimension.');
		end
		for i=1:q
			fprintf(fid,'k %d 0 %g\n',i,vlp.c(i)); 
		end
	end		

	% write rows  
	m1=max(size(vlp.a,1),size(vlp.b,1));
	if isempty(vlp.a)
		aa=-Inf*ones(m1,1);
	else
		aa=vlp.a;
	end
	if isempty(vlp.b)
		bb=Inf*ones(m1,1);
	else
		bb=vlp.b;
	end
	for i=1:m1
		if aa(i)<bb(i)
			ch=2*isfinite(aa(i))+isfinite(bb(i));
			switch ch
			case 0, fprintf(fid,'i %d f \n',i);
			case 1, fprintf(fid,'i %d u %g\n',i,bb(i));
			case 2, fprintf(fid,'i %d l %g\n',i,aa(i));
			case 3, fprintf(fid,'i %d d %g %g\n',i,aa(i),bb(i));
			end
		elseif aa(i)==bb(i) && isfinite(aa(i))
			fprintf(fid,'i %d s %g\n',i,aa(i));
		else
			error('Invalid constraints: a(%d)=%g, b(%d)=%g',i,aa(i),i,bb(i));
		end
	end

	% write cols
	if isempty(vlp.l)
		llb=-Inf*ones(n,1);
	else
		llb=vlp.l;
	end
	if isempty(vlp.s)
		uub=Inf*ones(n,1);
	else
		uub=vlp.s;
	end
	for j=1:n
		if llb(j)<uub(j)
			ch=2*isfinite(llb(j))+isfinite(uub(j));
			switch ch
			case 0, fprintf(fid,'j %d f\n',j); 
			case 1, fprintf(fid,'j %d u %g\n',j,uub(j));
			case 2, fprintf(fid,'j %d l %g\n',j,llb(j)); 
			case 3, fprintf(fid,'j %d d %g %g\n',j,llb(j),uub(j)); 
			end
		elseif llb(j)==uub(j) && isfinite(llb(j))
			fprintf(fid,'j %d s %g\n',j,llb(j));
		else
			error('Invalid constraints: l(%d)=%g, s(%d)=%g',j,llb(j),j,uub(j));
		end
	end 

	% end of file
	fprintf(fid,'e ');
	fclose(fid);
end