%--------------------------------------------------%
%          Compute first nm POD modes for U        %
%                                                  %
% Output: Mi(n1,n2,mode) -> spatial modes          %
%         Ci(mode,time)  -> temporal coefficients  %
%         Ei(mode)       -> associated energy      %
%                                                  %
%--------------------------------------------------%
function [Mu,Cu,Eu] = compute_POD(U,nm)

% get dimensions
n1     = size(U,1);
n2     = size(U,2);
ntimes = size(U,3);

% create matrix will all fluctuating velocities for each snapshot in a column 
Uall = reshape(U,n1*n2,ntimes);

% POD analysis
% Autocovariance matrix (ntimes x ntimes)
R = Uall'*Uall; 
% Compute eV, eigenvectors; D eigenvalues in diagonal matrix 
[eV,D] = eig(R);

% sort eigenvalues in ascending order. I is sorted index vector
[L,I]  = sort (diag( D));
eValue = zeros(size( L));
eVec   = zeros(size(eV)); % columns are eigenvectors. Eigenvector i: eVec(.,i)
for i=1:length(D)
    % Eigenvalues sorted in descending order
    eValue(length(D)+1-i) = L(i);       
    % Eigenvectors sorted in the same order
    eVec(:,length(D)+1-i) = eV(:,I(i)); 
end

% last eigenvalue should be zero 
eValue(length(eValue)) = 0;
 
% relative energy associated with mode m
Eu = eValue/sum(eValue); 

% calculate the first nm POD modes
Mu = zeros(n1*n2,nm);
for i=1:nm
    tmp     = Uall*eVec(:,i); % find mode. Sizes: (n1*n2,n)*(n,1)
    Mu(:,i) = tmp/norm(tmp);  % normalize mode 
end

% calculate POD coefficients
Cu = zeros(nm,ntimes);
for t=1:ntimes
    Cu(:,t) = Mu'*Uall(:,t);
end

% reshape modes
Mu = reshape(Mu,[n1 n2 nm]);

end
