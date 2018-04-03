R.<t> = PolynomialRing(QQ)

# the usual Verlinde formula
# = (1) in Zagier, Elementary aspects of the Verlinde formula and of the Harder-Narasimhan-Atiyah-Bott formula
def verlinde_sine(g, k):
  return round((k + 1)^(g - 1) * sum((-1)^(j-1) / (sin(j * pi / (2*k + 2))^(2*g-2)) for j in range(1, 2*k+2)))

# the determinantal Verlinde formula
# = (13) in Zagier, Elementary aspects of the Verlinde formula and of the Harder-Narasimhan-Atiyah-Bott formula
def verlinde_det(g, k):
  # we replace k by k+1 as (13) expresses D_-(g,2k), we need D_-(g,2k+2)
  # we take care that r=0,...,g-1 and s=1,...,g
  return 2^g * det(matrix(g, g, lambda r, s: 1 if r == 0 else (k+1+r)^(2*s+2) - (k+1-r)^(2*s+2))) / prod([factorial(2*i) for i in range(1, g+1)])

# Hilbert polynomial for M_C(2,L) with respect to Theta (not -K_X)
def Hilbert(g):
  U = sum([verlinde_det(g, k) * t^k for k in range(3*g + 1)]) * (1-t)^(3*g-2)
  return sum([U.coefficients(sparse=False)[k] * binomial(t - k + 3*(g-1), 3*(g-1)) for k in range(3*g)])


# checking (CS) for M_C(2,L)
for g in range(2, 20):
  print "(CS) for g=%d: %s" % (g, max([root[0].real() for root in Hilbert(g).roots(CC)]) <= 0)
print ""

# checking (NCS) for M_C(2,L)
for g in range(2, 20):
  print "(NCS) for g=%d: %s" % (g, max([root[0].real() for root in Hilbert(g).roots(CC)]) <= -2 / (3*g - 2))
print ""


# checking (CS) for hyperplane section of M_C(2,L)
def HilbertFano1(g): return Hilbert(g) - Hilbert(g).subs(t=t-1)
for g in range(2,20):
  print "(CS) for hyperplane section of M_C(2,L), g=%d: %s" % (g, max([root[0].real() for root in HilbertFano1(g).roots(CC)]) <= 0)
print ""

# checking (CS) for double cover of M_C(2,L) branched in Theta
def HilbertFano2(g): return Hilbert(g) + Hilbert(g).subs(t=t+1)
for g in range(2,20):
  print "(CS) for double cover of M_C(2,L) branched in Theta, g=%d: %s" % (g, max([root[0].real() for root in HilbertFano2(g).roots(CC)]) <= 0)
print ""


# checking (CL) for anticanonical section of M_C(2,L)
def HilbertCY1(g): return Hilbert(g) - Hilbert(g).subs(t=t-2)
for g in range(2,20):
  print "(CL) for anticanonical section of M_C(2,L), g=%d: %s" % (g, max([root[0].real() for root in HilbertCY1(g).roots(CC)]) == 0)
print ""

# checking (CL) for linear section of M_C(2,L) by pencil of hyperplanes
def HilbertCY2(g): return Hilbert(g) - 2*Hilbert(g).subs(t=t-1) + Hilbert(g).subs(t=t-2)
for g in range(2,20):
  print "(CL) for linear section of M_C(2,L) by pencil of hyperplanes, g=%d: %s" % (g, max([root[0].real() for root in HilbertCY2(g).roots(CC)]) == 0)
print ""

# checking (CL) for double cover of M_C(2,L) branched in -2K
def HilbertCY3(g): return Hilbert(g) + Hilbert(g).subs(t=t-2)
for g in range(2,20):
  print "(CL) for double cover of M_C(2,L) branched in -2K, g=%d: %s" % (g, max([root[0].real() for root in HilbertCY3(g).roots(CC)]) == 0)
print ""

# checking (CL) for cubic section of O(1)-cone of M_C(2,L)
def HilbertCY4(g): return Hilbert(g) + Hilbert(g).subs(t=t-1) + Hilbert(g).subs(t=t-2)
for g in range(2,20):
  print "(CL) for cubic section of O(1)-cone of M_C(2,L), g=%d: %s" % (g, max([root[0].real() for root in HilbertCY4(g).roots(CC)]) == 0)
print ""

# checking (CL) for double quartic section of double O(1)-cone of M_C(2,L)
def HilbertCY5(g): return Hilbert(g) + 2*Hilbert(g).subs(t=t-1) + Hilbert(g).subs(t=t-2)
for g in range(2,20):
  print "(CL) for double quartic section of double O(1)-cone of M_C(2,L), g=%d: %s" % (g, max([root[0].real() for root in HilbertCY5(g).roots(CC)]) == 0)
print ""

# checking (CL) for smoothing of a linear section of a join of M_C(2,L) with an elliptic curve of degree 1
def HilbertCY6(g): return Hilbert(g) - Hilbert(g).subs(t=t-1) + Hilbert(g).subs(t=t-2)
for g in range(2,20):
  print "(CL) for smoothing of a linear section of a join of M_C(2,L) with an elliptic curve of degree 1, g=%d: %s" % (g, max([root[0].real() for root in HilbertCY6(g).roots(CC)]) == 0)
