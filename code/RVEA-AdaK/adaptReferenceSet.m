function Wadapt = adaptReferenceSet(A, S, Wvalid, W, N, zmin, zmax, scale)
% Returns an adapted reference set of weight vectors
    zmax_arch = max(A, [],1);
    r = zmax_arch-zmin;
    epsilon = 1e-3;
    if (size(A,1) == N && all(r > epsilon))
        denom = zmax-zmin;
        denom(denom == 0) = 1e-12;
        V = (A(size(S,1)+1:end,:)-zmin)./denom;
        denom = sum(V, 1);
        denom(denom == 0) = 1e-12;
        Vweight = V./denom;
        Wadapt = [Wvalid;Vweight];
    else
        Wadapt = W;
    end
    if scale
        Wscale = Wadapt*(zmax-zmin);
        denom = sum(Wscale,1);
        denom(denom == 0) = 1e-12;
        Wadapt = Wscale./denom;
    end
end