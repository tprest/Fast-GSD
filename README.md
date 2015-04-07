Faster Algorithms for Gram-Schmidt Orthogonalization in NTRU lattices
===========


What is this library?
=====================
This software is a proof-of-concept Sage implementation of some of the algorithms described in [LP15]. The article can be found on http://www.di.ens.fr/~lyubash/ and http://www.di.ens.fr/~prest/.
More particularly, we focus on the Gram-Schmidt Orthogonalization (GSO) and Gram-Schmidt Decomposition (GSD) of NTRU lattices.
Since we use some clever ideas of [GHN06] to speed up our algorithms by taking advantage of the symplectic structure of NTRU lattices, we also implemented the GSD algorithms from [GHN06].


What does it contain?
=====================
The free open-source software Sage (aka SageMath) is necessary to run the algorithms implemented here.

The file NTRU.sage contains the tools to generate a public and a secret basis for a NTRU lattice.

The file GS.sage contains:
- Generic algorithms for performing the GSO/GSD over a generic field (Stable_GS) or the ring of integers (Standard_GS)
- Algorithms from [GHN06] for performing the GSD of generic bases (Dual_GS) or q-symplectic bases (Symplectic_GS) over the ring of integers

The file Isometric_GS.sage contains specialized algorithms from [LP15] for performing the GSO+GSD of NTRU lattices over generic field (Generic_NTRU_GSD) or the ring of integers (Integer_NTRU_GSD, Optimized_Integer_NTRU_GSD).


How do I use it?
================
On a Sage terminal, these commands allow to test the correctness of the GSO/GSD algorithms:
```
%attach NTRU.sage
%attach GS.sage
%attach Isometric_GS.sage

N=16        #Take N to be any power of 2
q=1024      #Take q to be any integer you want

(A,B) = Keygen(N,q);    #Generating a NTRU secret and public bases

#Testing the correctness of the rational GSD algorithms
print (Generic_NTRU_GSD(A,QQ,q) == A.gram_schmidt());

#Testing the correctness of the integer GSD algorithms
print (Standard_GS(A) == Dual_GS(A) == Symplectic_GS(A,q) == Integer_NTRU_GSD(A,q) == Optimized_Integer_NTRU_GSD(A,q) );
```


Possible improvements
=====================
If you have any suggestion or are interested on improving, adding new functionalities or porting the code, you can join me on my email address, which can be found on my website: http://www.di.ens.fr/~prest/


References
==========
[LP15] - Vadim Lyubashevsky and Thomas Prest: Quadratic Time, Linear Space Algorithms for Gram-Schmidt Orthogonalization and Gaussian Sampling in Structured Lattices

[GHN06] - Nicolas Gama, Nick Howgrave-Graham and Phong Q. Nguyen: Symplectic Lattice Reduction and NTRU
