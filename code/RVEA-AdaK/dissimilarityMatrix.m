function Diss = dissimilarityMatrix(A, flag)
%UNTITLED6 此处显示有关此函数的摘要
%   此处显示详细说明
    if(flag == 1)
        % m = np.shape(A)[1] 
        m=size(A,2);
        Diss = dissimilarityMatrixRSE(A, m-1);
    elseif(flag == 2)
        Diss = dissimilarityMatrixGAE(A, 512);
    elseif(flag == 3)
        Diss = dissimilarityMatrixCOU(A);
    elseif (flag == 4)
        Diss = dissimilarityMatrixPT(A, 5, 3, 0.02);
    elseif (flag == 5)
        Diss = dissimilarityMatrixMPT(A, 1, 25);
    elseif(flag == 6)
        Diss = dissimilarityMatrixKRA(A, 5, 3, 0.02);
    end
end

function d = dissimilarityMatrixRSE(A, s)
% Returns dissimilarity matrix using Riesz s-energy
    d = pdist(A);
    denom = d.^s;
    denom(denom == 0) = 1e-12;
    d = 1./denom;
    d=squareform(d);
end

function d = dissimilarityMatrixGAE(A, alpha)
% Returns dissimilarity matrix using Gaussian alpha-energy"""
    d = pdist(A);
    d = exp(-alpha.*(d.^2));
    d=squareform(d);
end

function d =dissimilarityMatrixCOU(A)
% Returns dissimilarity matrix using Coulomb's law"""
    k = 1/(4*pi*8.854e-12);
    norm = sum(A.^2, 2).^(1/2);
    % V = np.outer(norm, norm)
    V = norm*norm';
    V(logical(eye(size(V,1))))=0;
    % np.fill_diagonal(V, 0)
    % v = squareform(V);
    v=V;
    d = squareform(pdist(A, 'euclidean'));
    denom = d.^2;
    denom(denom == 0) = 1e-12;
    d = k*v./denom;
    % d=squareform(d);
end
function d= dissimilarityMatrixPT(A, V1, V2, alpha)
% Returns dissimilarity matrix using Pösch-Teller Potential"""
    d = pdist(A);
    denom1 = sin(alpha*d).^2;
    denom1(denom1 == 0) = 1e-12;
    denom2 = cos(alpha*d).^2;
    denom2(denom2 == 0) = 1e-12;
    d = V1/denom1+V2/denom2;
    d=squareform(d);
end

function d = dissimilarityMatrixMPT(A, D, alpha)
% Returns dissimilarity matrix using Modified Pösch-Teller Potential"""
    d = pdist(A);
    % d = -D/(cosh(alpha*d)**2)
    d=-D./(cosh(alpha*d).^(2));
    d=squareform(d);
end

function d= dissimilarityMatrixKRA(A, V1, V2, alpha)
% Returns dissimilarity matrix using Kratzer Potential"""
    d = pdist(A);
    denom = d;
    denom(denom == 0) = 1e-12;
    d = V1*(((d-(1/alpha))./denom).^2)+V2;
    d=squareform(d);
end