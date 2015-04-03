attach NTRU.sage


#The classical Gran-Schmidt algorithm (in truth a very common variant, known for its stability)
#This algorithm performs operations over the field F. One can typically take F to be QQ or RR.
def Stable_GS(A,F):
    N = A.nrows();
    M = A.ncols();
    B = matrix(F,N,M);
    v = vector(F,N);
    for i in range(N):
        B[i] = A[i];
        for j in range(i):
            B[i] = B[i] - (B[i]*B[j])/(v[j])*B[j];
        v[i] = B[i]*B[i];
    return B;



#Algorithm performing exact Gram-Schmidt Decomposition over generic lattices
#Transcribed from Algorithm 1 of [GHN06]
#This algorithm uses only integer arithmetic
def Standard_GS(B):
    G0 = B*transpose(B);
    n = B.nrows()//2;
    v = matrix(ZZ,1,2*n);
    w = matrix(ZZ, 1, 1);
    G = block_matrix([ [w, v], [transpose(v), G0] ]);
    U = vector(ZZ,2*n+1)
    Lamb = matrix(ZZ,2*n+1);

    Lamb[0,0] = 1;         #Sentinel value
    for i in [1..2*n]:
        Lamb[i,0] = 1;
        Lamb[i,1] = G[i,1];
        for j in [2..i]:
            S = Lamb[i,1]*Lamb[j,1];
            for k in [2..j-1]:
                S = (Lamb[k,k]*S + Lamb[j,k]*Lamb[i,k])/Lamb[k-1,k-1];
            Lamb[i,j] = G[i,j]*Lamb[j-1,j-1] - S;
    return Lamb.submatrix(1,1);



#Algorithm performing faster exact Gram-Schmidt Decomposition over symplectic lattices
#Transcribed from Algorithm 2 of [GHN06]
#This algorithm uses only integer arithmetic
def Dual_GS(B):
    G = B*transpose(B);
    n = B.nrows()//2;
    U = vector(ZZ,2*n+1)
    Lamb = zero_matrix(ZZ,2*n+1);

    for i in [0..2*n]:
        Lamb[i,0] = 1;         #Sentinel values

    for i in [1..2*n]:
        U[i] = Lamb[i-1,i-1];
        U[i-1] = -Lamb[i,i-1];
        li = [1..i-2];
        for j in reversed(li):
            U[j] = -sum(Lamb[k,j]*U[k] for k in [j+1..i])/Lamb[j,j];
        for j in [i..2*n]:
            Lamb[j,i] = sum(G[j-1,k-1]*U[k] for k in [1..i]);
    return Lamb.submatrix(1,1);



#Algorithm performing faster exact Gram-Schmidt Decomposition over symplectic lattices
#Transcribed from Algorithm 3 of [GHN06]
#This algorithm uses only integer arithmetic
def Symplectic_GS(B,q):
    G = B*transpose(B);
    n = B.nrows()//2;
    qq = [q^(2*i) for i in range(n+1)];
    U = vector(ZZ,2*n+1)
    Lamb = zero_matrix(ZZ,2*n+1);

    for i in [0..2*n]:
        Lamb[i,0] = 1;         #Sentinel values

    for i in [1..n]:
        U[i] = Lamb[i-1,i-1];
        U[i-1] = -Lamb[i,i-1];
        li = [1..i-2];
        for j in reversed(li):
            U[j] = -sum(Lamb[k,j]*U[k] for k in [j+1..i])/Lamb[j,j];
        for j in [1..i]:
            Lamb[2*n+1-j,2*n+1-i] = qq[n+1-i]*U[j];
        for j in [i..2*n]:
            Lamb[j,i] = sum(G[j-1,k-1]*U[k] for k in [1..i]);
    return Lamb.submatrix(1,1);
