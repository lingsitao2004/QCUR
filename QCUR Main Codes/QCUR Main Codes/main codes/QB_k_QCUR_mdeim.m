%time and error of QB_k_QCUR using deim and mdeim methods

% Reference: Efficient quaternion CUR decomposition based on discrete empirical interpolation method
% Sitao Ling and Zhehan Hu 2024

m=1024;n=1024;
A=qrandn(m,n);
relerr=1e-6; b=10;P=1; %opt=2;
t_mdeim=[];rel=[];t_deim=[];rel2=[];

for k=20:20:200

    
    Q = zeros(m, 0);
    B = zeros(0, n);
    r = 1;
    while r < k
        Omg = qrandn(n, b);
        Y = A * Omg - (Q * (B * Omg));
        [Qi,~] = QMGS_thinQR(Y);
        
        for j = 1:P        % power scheme
            [Qi,~] = QMGS_thinQR(A'*Qi - B'*(Q'*Qi)); % can skip orthonormalization for small b. 
            [Qi,~] = QMGS_thinQR(A*Qi - Q*(B*Qi));
        end
        
        if r>1,            % can skip the first re-orthogonalization
            [Qi,~] = QMGS_thinQR(Qi - Q * (Q' * Qi));
        end
        Bi= Qi'*A;  % another choice is Bi = Qi' * A - Qi' * Q * B;
        
        Q = [Q, Qi];
        B = [B; Bi];
        
        r = r + b;
    end
[U,S,V,~,~] = lansvdQ_restart([part(B,1) part(B,3) part(B,2) part(B,4)],k);
U=timesQ_vector([part(Q,1) part(Q,3) part(Q,2) part(Q,4)],U);
U=quaternion(U(:,1:k),U(:,2*k+1:3*k),U(:,k+1:2*k),U(:,3*k+1:4*k));
V=quaternion(V(:,1:k),V(:,2*k+1:3*k),V(:,k+1:2*k),V(:,3*k+1:4*k));

%% mdeim

irow=zeros(1,k);
icol=zeros(1,k);
tic
for j = 1:k
  [~, irow(j)] = max(abs(U(:,j)));
  [~, icol(j)] = max(abs(V(:,j)));
  if j<k 
      if j==1
      zz1=pinvq(U(irow(1:j),1:j));
      zz2=pinvq(V(icol(1:j),1:j));
      U(:,j+1) = U(:,j+1) - U(:,1:j) * zz1*U(irow(1:j),j+1);
      V(:,j+1) = V(:,j+1) - V(:,1:j) * zz2*V(icol(1:j),j+1);
      else
      zz1=sher_inv(zz1,U(irow(1:j-1),j),U(irow(j),1:j-1),U(irow(j),j));
      zz2=sher_inv(zz2,V(icol(1:j-1),j),V(icol(j),1:j-1),V(icol(j),j));

      U(:,j+1) = U(:,j+1) - U(:,1:j) * zz1*U(irow(1:j),j+1);
      V(:,j+1) = V(:,j+1) - V(:,1:j) * zz2*V(icol(1:j),j+1);
      end
  end
end
 t1=toc;
 t_mdeim=[t_mdeim,t1]
 C=A(:,icol);
 R=A(irow,:);
 CP=pinv(C'*C);CP=CP*C';
 RP=pinv(R*R');RP=R'*RP;
 M=A*RP;
 M=CP*M;

rel=[rel,norm(C*M*R-A,'fro')/norm(A,'fro')]

%% deim
  tic
  for j = 1:k
    [~, irow(j)] = max(abs(U(:,j)));
    [~, icol(j)] = max(abs(V(:,j)));
    if j<k
     U(:,j+1) = U(:,j+1) - U(:,1:j) * pinvq(U(irow(1:j),1:j))*U(irow(1:j),j+1);
     V(:,j+1) = V(:,j+1) - V(:,1:j) * pinvq(V(icol(1:j),1:j))*V(icol(1:j),j+1);
    end
  end
  t2=toc;
  t_deim=[t_deim,t2]
  C=A(:,icol);
  R=A(irow,:);
  CP=pinv(C'*C);CP=CP*C';
  RP=pinv(R*R');RP=R'*RP;
  M=A*RP;
  M=CP*M;

  rel2=[rel2,norm(C*M*R-A,'fro')/norm(A,'fro')]

end

