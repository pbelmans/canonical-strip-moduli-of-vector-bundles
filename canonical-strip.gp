\\ increase precision
\p 100

\\ MaxRealRoot(p) returns maximal real value of a root of polynomial p
MaxRealRoot(p) = vecmax(real(polroots(p)));

\\ the usual Verlinde formula
\\ = (1) in Zagier, Elementary aspects of the Verlinde formula and of the Harder-Narasimhan-Atiyah-Bott formula
\\ beware of rounding errors
verlinde_sine(g,k) = round((k+1)^(g-1)*sum(j=1,2*k+1,(-1)^(j-1)/sin(Pi*j/(2*k+2))^(2*g-2)));

\\ the determinantal Verlinde formula
\\ = (13) in Zagier, Elementary aspects of the Verlinde formula and of the Harder-Narasimhan-Atiyah-Bott formula
\\ we replace k by k+1 as (13) expresses D_-(g,2k), we need D_-(g,2k+2)
\\ we take care that r=0,...,g-1 and s=1,...,g
verlinde_det(g,k) = 2^g*matdet(matrix(g,g,r,s,if(r==1,1,(k+r)^(2*s)-(k-r+2)^(2*s))))/prod(i=1,g,(2*i)!,1)

\\ numerator of Hilbert--Poincare series for M
U(g) = truncate(sum(n=0,3*g, verlinde_det(g,n) * t^n, O(t^(3*g+1)))*(1-t)^(3*g-2));

\\ Hilbert polynomial for M_C(2,L) with respect to Theta (not -K_X)
Hilbert(g) = {UU=U(g);return(sum(k=0,poldegree(UU),polcoeff(UU,k)*binomial(n-k+3*(g-1),3*(g-1))))};

\\ checking (CS) for M_C(2,L)
for(g=2,20,print("(CS) for M_C(2,L), g="g,": ",MaxRealRoot(Hilbert(g)) <= 0));

\\ checking (NCS) for M_C(2,L)
for(g=2,20,print("(NCS) for M_C(2,L), g="g,": ",MaxRealRoot(Hilbert(g)) <= -2 / (3*g-2)));


\\ checking (CS) for a hyperplane section of M_C(2,L)
HilbertFano1(g) = {HH=Hilbert(g); return(HH-subst(HH,n,n-1))};
for(g=2,20,print("(CS) for hyperplane section of M_C(2,L), g="g,": ",MaxRealRoot(HilbertFano1(g)) <= 0));

\\ checking (CS) for double cover of M_C(2,L) branched in Theta
HilbertFano2(g) = {HH=Hilbert(g); return(HH+subst(HH,n,n-1))};
for(g=2,20,print("(CS) for double cover of M_C(2,L) branched in Theta, g="g,": ",MaxRealRoot(HilbertFano2(g)) <= 0));


\\ checking (CL) for anticanonical section of M_C(2,L)
HilbertCY1(g) = {HH=Hilbert(g); return(HH-subst(HH,n,n-2))};
for(g=2,20,print("(CL) for anticanonical section of M_C(2,L), g="g,": ",MaxRealRoot(HilbertCY1(g)) == 0));

\\ checking (CL) for linear section of M_C(2,L) by pencil of hyperplanes
HilbertCY2(g) = {HH=Hilbert(g); return(HH-2*subst(HH,n,n-1)+subst(HH,n,n-2))};
for(g=2,20,print("(CL) for linear section of M_C(2,L) by pencil of hyperplanes, g="g,": ",MaxRealRoot(HilbertCY2(g)) == 0));

\\ checking (CL) for double cover of M_C(2,L) branched in -2K
HilbertCY3(g) = {HH=Hilbert(g); return(HH+subst(HH,n,n-2))};
for(g=2,20,print("(CL) for double cover of M_C(2,L) branched in -2K, g="g,": ",MaxRealRoot(HilbertCY3(g)) == 0));

\\ checking (CL) for cubic section of O(1)-cone of M_C(2,L)
HilbertCY4(g) = {HH=Hilbert(g); return(HH+subst(HH,n,n-1)+subst(HH,n,n-2))};
for(g=2,20,print("(CL) for cubic section of O(1)-cone of M_C(2,L), g="g,": ",MaxRealRoot(HilbertCY4(g)) == 0));

\\ checking (CL) for double quartic section of double O(1)-cone of M_C(2,L)
HilbertCY5(g) = {HH=Hilbert(g); return(HH+2*subst(HH,n,n-1)+subst(HH,n,n-2))};
for(g=2,20,print("(CL) for double quartic section of double O(1)-cone of M_C(2,L), g="g,": ",MaxRealRoot(HilbertCY5(g)) == 0));

\\ checking (CL) for smoothing of a linear section of a join of M_C(2,L) with an elliptic curve of degree 1
HilbertCY6(g) = {HH=Hilbert(g); return(HH-subst(HH,n,n-1)+subst(HH,n,n-2))};
for(g=2,20,print("(CL) for smoothing of a linear section of a join of M_C(2,L) with an elliptic curve of degree 1, g="g,": ",MaxRealRoot(HilbertCY6(g)) == 0));
