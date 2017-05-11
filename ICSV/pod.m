function[phi, lmd, R] = pod(V)
%	Snapshot POD assumes each row contains all measurements taken at a given time
	[Nx, Nt] = size(V);
	wgt = zeros(Nx,Nt);
%	for i = 1:Nx
%		for j = 1:Nt
%			wgt(i,j) = i;
%			V(i,j) = V(i,j)*wgt(i,j);
%		end
%		V(i,:) = V(i,:) - mean(V(i,:));
%	end
%	V = V.*wgt;


	U = V';
%	1) Compute the spatial correlation matrix
	C = U*U';
%	2) Solve the eigenvalue problem
	[beta, lmd] = svd(C);
%	3) Calculate and normeralize the basis functions
	phi = U'*beta;
%	phi = phi./wgt;
	normer = 0;
	for j = 1:Nt
		for i = 1:Nx
			normer = normer + phi(i,j).^2;
		end
	end
	normer = sqrt(normer);
	phi = phi./normer;
%	4) Compute the temporal coefficients
	R = U*phi;
end
