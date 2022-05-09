function [Mu,Phi] = compute_DMD(X,Y,nm)
       
    [U,S,V] = svd(X);
    
    U = U(:,1:nm);
    S = S(1:nm,1:nm);
    V = conj(V);
    V = conj(V(:,1:nm));
    
    A_tilde = ((conj(U)'*Y)*V)*inv(S);
    
    [W,Mu] = eig(A_tilde);
    
    Phi = Y*V*inv(S)*W;
    
    %Phi = reshape(Phi,[46,33,nm]);
    %Phi = permute(Phi,[2,1,3]);

end