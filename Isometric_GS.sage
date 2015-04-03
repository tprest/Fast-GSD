attach NTRU.sage


#Rotate k times a NTRU lattice vector
#The rotation r is defined as follows: r(f,g) = (x*f mod (x^n+1), x*g mod (x^n+1))
def NTRU_Rotate(v,F,k):
#This algorithm performs operations over F. One can typically take F to be ZZ, QQ or RR.
    n = len(v)//2;
    w = vector(F, 2*n);
    for i in [0..n-1-k]:
        w[i+k] = v[i];
        w[n+i+k] = v[n+i];
    for i in range(k):
        w[i] = -v[n-k+i];
        w[n+i] = -v[2*n-k+i];
    return w;


#Algorithm performing Gram-Schmidt Decomposition over NTRU lattices using the principles of [GHN06] and [LP15].
#Optimized versions of Algorithms 3,4 of [LP15]. See also section 7 of [LP15].
#This algorithm performs operations over the field F. One can typically take F to be QQ or RR.
def Generic_NTRU_GSD(A,F,q):
    n = A.nrows()//2;
    assert(A.nrows()==A.ncols());
    B = matrix(F,2*n,2*n);
    C = vector(F,2*n);
    D = vector(F,2*n);
    CD = vector(F,2*n);
    M = matrix(F,2*n,2*n);
    X = matrix(F,2*n,2*n);
    Y = matrix(F,2*n,2*n);

    B[0] = A[0];
    v = A[0];
    C[0] = NTRU_Rotate(B[0],F,1)*v;
    D[0] = v*v;
    CD[0] = C[0]/D[0];

    for i in [1..n-1]:
        B[i] = NTRU_Rotate(B[i-1],F,1) - CD[i-1]*v;
        v = v - CD[i-1]*NTRU_Rotate(B[i-1],F,1);
        Y[n,i] = NTRU_Rotate(A[n],F,1)*v;    #This is computed now to be used later in the GSD.
        D[i] = D[i-1] - C[i-1]*CD[i-1];
        C[i] =  NTRU_Rotate(B[i],F,1)*A[0];
        CD[i] = C[i]/D[i];

    for i in range(n):
        k = q/D[i];

        for j in range(n):
            B[2*n-1-i,2*n-1-j] = k*B[i,j];
            B[2*n-1-i,j] = -k*B[i,2*n-1-j];
        D[2*n-1-i] = q*k;

    vv = B[n,1];
    for i in range(n-1):
        assert(vv!=0);
        CD[n+i] = (B[n+i,0]-B[n+i+1,1])/vv;
        vv = vv - CD[n+i]*B[n+i,0];


    for i in range(2*n):
        M[i,i] = 1;


    #We fill the upper-left quadrant:
    for i in range(n):
        X[i,i] = D[i];
    for i in [1..n-1]:
        X[i,0] = A[i]*B[0];
        Y[i,0] = NTRU_Rotate(A[i],F,1)*A[0];
        for j in [1..i-1]:
            X[i,j] = X[i-1,j-1] - CD[j-1]*Y[i-1,j-1];
            Y[i,j] = Y[i,j-1] - CD[j-1]*X[i,j-1];


    #We fill the lower-left quadrant (under its diagonal):
    for i in [n..2*n-1]:
        X[i,0] = A[i]*B[0];
        Y[i,0] = NTRU_Rotate(A[i],F,1)*A[0];
        for j in [1..i-n]:
            X[i,j] = X[i-1,j-1] - CD[j-1]*Y[i-1,j-1];
            Y[i,j] = Y[i,j-1] - CD[j-1]*X[i,j-1];


    #We fill the lower-left quadrant (over its diagonal):
    for j in range(n):
        X[n,j] = A[n]*B[j];
        for i in [n+1..j+n-1]:
            X[i,j] = X[i-1,j-1] - CD[j-1]*Y[i-1,j-1];
            Y[i,j] = Y[i,j-1] - CD[j-1]*X[i,j-1];


    #We fill the lower-right quadrant:
    for i in [n..2*n-1]:
        X[i,n] = A[i]*B[n];
        Y[i,n] = NTRU_Rotate(A[i],F,1)*B[n];
        for j in [n+1..i]:
            X[i,j] = X[i-1,j-1] - CD[j-1]*Y[i-1,j-1];
            Y[i,j] = Y[i,j-1] - CD[j-1]*X[i,j-1];


    for i in range(2*n):
        for j in [0..i-1]:
            M[i,j] = X[i,j]/D[j];


    return (B,M);



#Algorithm performing Gram-Schmidt Decomposition over NTRU lattices using the principles of [GHN06] and [LP15].
#Optimized versions of Algorithms 6,7 of [LP15]. See also section 7 of [LP15].
#This algorithm performs only exact operations over the ring of integers ZZ
def Integer_NTRU_GSD(A,q):
    n = A.nrows()//2;

    b = matrix(ZZ,2*n,2*n);
    c = vector(ZZ,2*n);
    d = vector(ZZ,2*n+1);
    xx = matrix(ZZ,2*n,2*n);
    y = matrix(ZZ,2*n,2*n);

    d[2*n] = 1;     #This is the sentinel value for 'd[-1]'
    b[0] = A[0];
    v = A[0];
    c[0] = NTRU_Rotate(A[0],ZZ,1)*v;
    d[0] = v*v;

    for i in [1..n-1]:
        b[i] = (d[i-1]*NTRU_Rotate(b[i-1],ZZ,1) - c[i-1]*v)/d[i-2];
        v = (d[i-1]*v - c[i-1]*NTRU_Rotate(b[i-1],ZZ,1))/d[i-2]; 
        y[n,i] = NTRU_Rotate(A[n],ZZ,1)*v;                   #This is computed now to be used later in the GSD.
        d[i] = (d[i-1]*d[i-1] - c[i-1]*c[i-1])/d[i-2];
        c[i] = A[0]*NTRU_Rotate(b[i],ZZ,1);


    #Using symplecticity (see [GHN06]) to compute the second half of the GSO
    q2 = 1;
    for i in range(n):
        q2 *= (q**2);
        d[n+i] = q2*d[n-2-i];


    #Using symplecticity (see [GHN06]) to compute the second half of the GSO
    for i in range(n):
        k = q**(2*n-2*i-1);
        for j in range(n):
            b[2*n-1-i,2*n-1-j] = k*b[i,j];
            b[2*n-1-i,j] = -k*b[i,2*n-1-j];


    #Can be computed faster (see Generic_NTRU_GSD)
    for i in range(n-1):
        c[n+i] = (b[N]*NTRU_Rotate(b[n+i],ZZ,1))/d[n-1];


    #We fill the diagonal of xx:
    for i in range(2*n):
        xx[i,i] = d[i];


    #We fill the upper-left quadrant of xx:
    for i in [1..N-1]:
        xx[i,0] = A[i]*b[0];
        y[i,0] = NTRU_Rotate(A[i],ZZ,1)*A[0];
        for j in [1..i-1]:
            xx[i,j] = (d[j-1]*xx[i-1,j-1] - c[j-1]*y[i-1,j-1])/d[j-2];
            y[i,j] = (d[j-1]*y[i  ,j-1] - c[j-1]*xx[i  ,j-1])/d[j-2];


    #We fill the lower-left quadrant of xx (under its diagonal):
    for i in [n..2*n-1]:
        xx[i,0] = A[i]*b[0];
        y[i,0] = NTRU_Rotate(A[i],ZZ,1)*A[0];
        for j in [1..i-n]:
            xx[i,j] = (d[j-1]*xx[i-1,j-1] - c[j-1]*y[i-1,j-1])/d[j-2];
            y[i,j] = (d[j-1]*y[i  ,j-1] - c[j-1]*xx[i  ,j-1])/d[j-2];


    #We fill the lower-left quadrant of xx (over its diagonal):
    for j in range(n):
        xx[n,j] = A[n]*b[j];
        for i in [n+1..j+n-1]:
            xx[i,j] = (d[j-1]*xx[i-1,j-1] - c[j-1]*y[i-1,j-1])/d[j-2];
            y[i,j] = (d[j-1]*y[i  ,j-1] - c[j-1]*xx[i  ,j-1])/d[j-2];


    #We fill the lower-right quadrant of xx:
    for i in [n..2*n-1]:
        xx[i,n] = A[i]*b[n];
        y[i,n] = NTRU_Rotate(A[i],ZZ,1)*b[n];
        for j in [n+1..i-1]:
            xx[i,j] = (d[j-1]*xx[i-1,j-1] - c[j-1]*y[i-1,j-1])/d[j-2];
            y[i,j] = (d[j-1]*y[i  ,j-1] - c[j-1]*xx[i  ,j-1])/d[j-2];


    return xx;



#Allows to compute multiple dot products at the same time when the second term is constant and the first term is NTRU_Rotate(a,i) for i in range(n)
#Useful for algorithm Optimized_Integer_NTRU_GSD
def Multiple_Dot_Products(a,b,F):
    assert(len(a)==len(b));
    n = len(a)//2;
    phi = x^n + 1;
    a1 = sum(a[i  ]*(x**i) for i in range(n));
    a2 = sum(a[i+n]*(x**i) for i in range(n));
    b1 = sum(b[i  ]*(x**i) for i in range(n));
    b2 = sum(b[i+n]*(x**i) for i in range(n));
    a1 = Reverse(a1,n);
    a2 = Reverse(a2,n);
    a1 = a1.change_ring(F);
    a2 = a2.change_ring(F);
    b1 = b1.change_ring(F);
    b2 = b2.change_ring(F);
    return (a1*b1 + a2*b2)%phi;


#Algorithm performing Gram-Schmidt Decomposition over NTRU lattices using the principles of [GHN06] and [LP15].
#Optimized versions of Algorithms 6,7 of [LP15]. See also section 7 of [LP15].
#This algorithm performs only exact operations over the ring of integers ZZ
#This algorithm introduces further optimizations based on the fact that given certain conditions, you can compute n dot products at essentially the cost of one dot product
def Optimized_Integer_NTRU_GSD(A,q):
    n = A.nrows()//2;

    b = matrix(ZZ,2*n,2*n);
    c = vector(ZZ,2*n);
    d = vector(ZZ,2*n+1);
    xx = matrix(ZZ,2*n,2*n);
    y = matrix(ZZ,2*n,2*n);


    d[2*n] = 1;     #This is the sentinel value for 'd[-1]'
    b[0] = A[0];
    v = A[0];
    c[0] = NTRU_Rotate(A[0],ZZ,1)*v;
    d[0] = v*v;


    for i in [1..n-1]:
        b[i] = (d[i-1]*NTRU_Rotate(b[i-1],ZZ,1) - c[i-1]*v)/d[i-2];
        v = (d[i-1]*v - c[i-1]*NTRU_Rotate(b[i-1],ZZ,1))/d[i-2]; 
        y[n,i] = NTRU_Rotate(A[n],ZZ,1)*v;                   #This is computed now to be used later in the GSD.
        d[i] = (d[i-1]*d[i-1] - c[i-1]*c[i-1])/d[i-2];
        c[i] = A[0]*NTRU_Rotate(b[i],ZZ,1);


    #Using symplecticity (see [GHN06]) to compute the second half of d
    q2 = 1;
    for i in range(n):
        q2 *= (q**2);
        d[n+i] = q2*d[n-2-i];


    #Using symplecticity (see [GHN06]) to compute the second half of the GSO
    for i in range(n):
        k = q**(2*n-2*i-1);
        for j in range(n):
            b[2*n-1-i,2*n-1-j] = k*b[i,j];
            b[2*n-1-i,j] = -k*b[i,2*n-1-j];


    #Can be computed faster (see Generic_NTRU_GSD)
    for i in range(n-1):
        c[n+i] = (b[n]*NTRU_Rotate(b[n+i],ZZ,1))/d[n-1];


    #We fill the diagonal of xx:
    for i in range(2*n):
        xx[i,i] = d[i];


    #We fill the upper-left quadrant of xx:
    MDP1 = Multiple_Dot_Products(A[0],A[0],ZZ);
    MDP2 = Multiple_Dot_Products(NTRU_Rotate(A[0],ZZ,1),A[0],ZZ);
    for i in [1..N-1]:
        xx[i,0] = MDP1[i];
        y[i,0] = MDP2[i];
        for j in [1..i-1]:
            xx[i,j] = (d[j-1]*xx[i-1,j-1] - c[j-1]*y[i-1,j-1])/d[j-2];
            y[i,j]  = (d[j-1]*y[i  ,j-1] - c[j-1]*xx[i  ,j-1])/d[j-2];


    #We fill the lower-left quadrant of xx (under its diagonal):
    MDP1 = Multiple_Dot_Products(A[n],A[0],ZZ);
    MDP2 = Multiple_Dot_Products(NTRU_Rotate(A[n],ZZ,1),A[0],ZZ);
    for i in [n..2*n-1]:
        xx[i,0] = MDP1[i-n];
        y[i,0] = MDP2[i-n];
        for j in [1..i-n]:
            xx[i,j] = (d[j-1]*xx[i-1,j-1] - c[j-1]*y[i-1,j-1])/d[j-2];
            y[i,j]  = (d[j-1]*y[i  ,j-1] - c[j-1]*xx[i  ,j-1])/d[j-2];


    #We fill the lower-left quadrant of xx (over its diagonal):
    for j in range(n):
        xx[n,j] = A[n]*b[j];
        for i in [n+1..j+n-1]:
            xx[i,j] = (d[j-1]*xx[i-1,j-1] - c[j-1]*y[i-1,j-1])/d[j-2];
            y[i,j]  = (d[j-1]*y[i  ,j-1] - c[j-1]*xx[i  ,j-1])/d[j-2];


    MDP1 = Multiple_Dot_Products(A[n],b[n],ZZ);
    MDP2 = Multiple_Dot_Products(NTRU_Rotate(A[n],ZZ,1),b[n],ZZ);
    #We fill the lower-right quadrant of xx:
    for i in [n..2*n-1]:
        xx[i,n] = MDP1[i-n];
        y[i,n] = MDP2[i-n];
        for j in [n+1..i-1]:
            xx[i,j] = (d[j-1]*xx[i-1,j-1] - c[j-1]*y[i-1,j-1])/d[j-2];
            y[i,j]  = (d[j-1]*y[i  ,j-1] - c[j-1]*xx[i  ,j-1])/d[j-2];


    return xx;
