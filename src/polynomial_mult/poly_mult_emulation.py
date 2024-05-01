# Python script with a few utilities to experiment with polynomial multiplication
import random

class Ring:
    """ Extended Ring structure """
    def __init__(self, modulo, degree, rootOfUnity, invRootOfUnity, invDegree):
        self.modulo = modulo
        self.degree = self.elt(degree)
        self.invDegree = self.elt(invDegree)
        self.rootOfUnity = self.elt(rootOfUnity)
        self.invRootOfUnity = self.elt(invRootOfUnity)
        if not invRootOfUnity is None: assert (rootOfUnity * invRootOfUnity) % modulo == 1
        assert (self.rootOfUnity ** degree).value % modulo == 1
        if not invDegree is None: assert (degree * invDegree) % modulo == 1

    def elt(self, value):
        """ build and encapsulate an element to the Ring @p self """
        if isinstance(value, RingElt):
            # assert value.ring == self
            return value
        if isinstance(value, int):
            return RingElt(value, self)
        if value is None:
            return None
        print(value.__class__)
        raise NotImplementedError

    def powerRing(self):
        """ sub-ring with half the degreee and the square of the root of unity"""
        return Ring(self.modulo, self.degree.value // 2, self.rootOfUnity ** 2, None, None)

    def __eq__(lhs, rhs) -> bool:
        """ testing Ring equality """
        # limiting equality check to modulo, degree and rootOfUnity
        return lhs.modulo == rhs.modulo and lhs.rootOfUnity == rhs.rootOfUnity and lhs.degree == rhs.degree

    @property
    def zero(self):
        """ build element zero of Ring @p self """
        return RingElt(0, self)

    def __repr__(self) -> str:
        return f"Ring Z/{self.modulo}.Z"


class RingElt:
    """ Ring element """
    def __init__(self, value, ring):
        self.value = value
        self.ring = ring

    def __add__(lhs, rhs):
        assert lhs.ring == rhs.ring
        return RingElt((lhs.value + rhs.value) % lhs.ring.modulo, lhs.ring)

    def __sub__(lhs, rhs):
        assert lhs.ring == rhs.ring
        return RingElt((lhs.value - rhs.value) % lhs.ring.modulo, lhs.ring)

    def __mul__(lhs, rhs):
        assert lhs.ring == rhs.ring
        return RingElt((lhs.value * rhs.value) % lhs.ring.modulo, lhs.ring)

    def __pow__(self, n):
        if n == 0:
            return RingElt(1, self.ring)
        if n == 1:
            return self
        if n > 1:
            return self * self**(n-1)
        raise NotImplementedError("power must be positive")

    def __neg__(self):
        return RingElt(-self.value, self.ring)

    def __repr__(self) -> str:
        return f"{self.value}"

    def __eq__(self, value: object) -> bool:
        if isinstance(value, int):
            return self.value == value
        elif isinstance(value, RingElt):
            return self.value == value.value
        raise NotImplementedError
        

class Polynomial:
    """ Structure to encode a Polynomial object """
    def __init__(self, coeffs, eltRing: Ring):
        self.degree = len(coeffs) - 1
        self.coeffs = [eltRing.elt(v) for v in coeffs]
        self.eltRing = eltRing

    def __add__(lhs, rhs):
        assert lhs.eltRing == rhs.eltRing
        tailStart = min(lhs.degree, rhs.degree) + 1
        # at most one of rhs.coeffs[tailStart:] and lhs.coeffs[tailStart:] is non-empty
        coeffs = [l + r for (l, r) in zip(lhs.coeffs, rhs.coeffs)] + rhs.coeffs[tailStart:] + lhs.coeffs[tailStart:]
        return Polynomial(coeffs, lhs.eltRing)

    def __sub__(lhs, rhs):
        assert lhs.eltRing == rhs.eltRing
        if lhs.degree >= rhs.degree:
            tail = lhs.coeffs[rhs.degree + 1:]
        else:
            tail = [-v for v in rhs.coeffs[lhs.degree + 1:]]
        coeffs = [l - r for (l, r) in zip(lhs.coeffs, rhs.coeffs)] + tail
        return Polynomial(coeffs, lhs.eltRing)


    def __mul__(lhs, rhs):
        assert lhs.eltRing == rhs.eltRing
        coeffs = [lhs.eltRing.zero] * (lhs.degree + rhs.degree + 1) # over provisioned if both degress are 0
        for li, lv in enumerate(lhs.coeffs):
            for ri, rv in enumerate(rhs.coeffs):
                coeffs[li + ri] += lv * rv
        return Polynomial(coeffs, lhs.eltRing)

    def __mod__(self, mod):
        assert(mod.coeffs[-1] == 1) # largest coefficients must be one
        def helperMod(lhs):
            if lhs.degree < mod.degree:
                return lhs
            else:
                scale = lhs.degree - mod.degree
                factor = Polynomial.constant(lhs.coeffs[-1], self.eltRing)
                remainder = lhs - factor * Polynomial.monomial(scale, self.eltRing) * mod
                assert(remainder.coeffs[-1] == 0)
                return helperMod(Polynomial(remainder.coeffs[:-1], self.eltRing))
        return helperMod(self)

    def __repr__(self):
        if self.isZero():
            return "0"
        cst = [] if self.coeffs[0] == 0 else [f"{self.coeffs[0]}"]
        degree1 = [] if self.degree < 1 else ([] if self.coeffs[1] == 0 else ["X" if self.coeffs[1] == 1 else f"{self.coeffs[1]}.X"])
        return " + ".join(cst + degree1 + [f"X^{i}" if v == 1 else f"{v}.X^{i}" for (i, v) in filter((lambda v: v[1] != 0), enumerate(self.coeffs[2:], 2))])

    def isZero(self):
        return all(v == 0 for v in self.coeffs)

    def extractEven(self):
        """ extract sub-polynomial with even-index coefficients """
        return Polynomial([self.coeffs[i] for i in range(0, self.degree + 1, 2)], self.eltRing.powerRing())

    def extractOdd(self):
        """ extract sub-polynomial with odd-index coefficients """
        return Polynomial([self.coeffs[i] for i in range(1, self.degree + 1, 2)], self.eltRing.powerRing())

    @staticmethod
    def monomial(degree, eltRing):
        return Polynomial([0] * degree + [1], eltRing)

    @staticmethod
    def constant(value, eltRing):
        return Polynomial([value], eltRing)

    @staticmethod
    def random(degree, eltRing):
        return Polynomial([random.randrange(eltRing.modulo) for d in range(degree + 1)], eltRing)

class NTTDomain:
    """ Structure to encode a value in the NTT domain """
    def __init__(self, coeffs, eltRing):
        self.degree = len(coeffs) - 1
        self.coeffs = [eltRing.elt(v) for v in coeffs]
        self.eltRing = eltRing

    def __add__(lhs, rhs):
        assert lhs.eltRing == rhs.eltRing
        assert len(lhs.coeffs) == len(rhs.coeffs)
        coeffs = [l + r for (l, r) in zip(lhs.coeffs, rhs.coeffs)]
        return Polynomial(coeffs, lhs.eltRing)

    def __sub__(lhs, rhs):
        assert lhs.eltRing == rhs.eltRing
        assert len(lhs.coeffs) == len(rhs.coeffs)
        coeffs = [l - r for (l, r) in zip(lhs.coeffs, rhs.coeffs)]
        return Polynomial(coeffs, lhs.eltRing)

    def __mul__(lhs, rhs):
        # assume the NTT domain allows element-wise multiplication
        assert lhs.eltRing == rhs.eltRing
        assert len(lhs.coeffs) == len(rhs.coeffs)
        coeffs = [l * r for (l, r) in zip(lhs.coeffs, rhs.coeffs)]
        return Polynomial(coeffs, lhs.eltRing)


def ntt_transform(poly: Polynomial) -> NTTDomain:
    """ Number Theoretic transform on @p poly """
    ntt_coeffs = []
    for k in range(len(poly.coeffs)):
        ntt_coeff = poly.eltRing.elt(0)
        for j in range(len(poly.coeffs)):
            ntt_coeff += poly.coeffs[j] * poly.eltRing.rootOfUnity ** (j * k)
        ntt_coeffs.append(ntt_coeff)
    return NTTDomain(ntt_coeffs, poly.eltRing)

def ntt_inv_transform(ntt: NTTDomain) -> Polynomial:
    """ Inverse Number Theoretic Transform """
    poly_coeffs = []
    for k in range(len(ntt.coeffs)):
        poly_coeff = ntt.eltRing.elt(0)
        for j in range(len(ntt.coeffs)):
            poly_coeff += ntt.coeffs[j] * ntt.eltRing.invRootOfUnity ** (j * k)
        poly_coeffs.append(ntt.eltRing.invDegree * poly_coeff)
    return Polynomial(poly_coeffs, ntt.eltRing)

def fast_ntt_transform(poly: Polynomial, inv: bool = False) -> NTTDomain:
    """ fast implementation of NTT (based on Cooler-Tukey FFT algorithm) """
    if poly.degree == 0:
        # print("fast NTT: level 1 reached")
        return NTTDomain(poly.coeffs[:1], poly.eltRing)
    if poly.degree % 2 == 0:
        # print(f"fast NTT: even degree {poly.degree}, fall back to slow NTT")
        # only works if the polynomial has an even number of coefficients
        # else fallback to slow method
        return ntt_transform(poly)
    else:
        # print("fast NTT: split and reconstruction")
        poly_even = poly.extractEven()
        poly_odd = poly.extractOdd()
        # recursive NTT
        ntt_even = fast_ntt_transform(poly_even)
        ntt_odd = fast_ntt_transform(poly_odd)
        # reconstruction
        ntt_coeffs = []
        half_degree = (poly.degree + 1) // 2
        for i in range(poly.degree + 1):        
            j = i % half_degree
            even_coeff = ntt_even.coeffs[j] 
            odd_coeff = ntt_odd.coeffs[j] 
            twiddle = (poly.eltRing.invRootOfUnity if inv else poly.eltRing.rootOfUnity) ** j
            if i >= half_degree:
                twiddle = -twiddle
            ntt_coeff = even_coeff + twiddle * odd_coeff
            ntt_coeffs.append(ntt_coeff)
        return NTTDomain(ntt_coeffs, poly.eltRing)


def fast_inv_ntt_transform(ntt: NTTDomain) -> Polynomial:
    """ Fast implementation of the inverse Number Theoretic Transform (NTT) """
    # re-using the Fast NTT implementation: transforming the input NTT representation
    # to a dummy polynomial so we can execute the fast NTT transform and then transforming back
    # the result representation in NTT domain to the output polynomial (while scaling by 1 / (degree+1))
    dummy_poly = Polynomial(ntt.coeffs, ntt.eltRing) 
    pre_poly_coeffs = fast_ntt_transform(dummy_poly, inv=True)
    return Polynomial([c * dummy_poly.eltRing.invDegree for c in pre_poly_coeffs.coeffs], ntt.eltRing)


if __name__ == "__main__":
    # DefaultRing = Ring(19, 3, 11, 7, 13)
    # declaring basis polynomial
    # mod = Polynomial.monomial(3, DefaultRing) + Polynomial.constant(-1, DefaultRing)
    DefaultRing = Ring(3329, 4, (17 ** 64) % 3329, 1600, 2497)
    mod = Polynomial.monomial(4, DefaultRing) + Polynomial.constant(-1, DefaultRing)

    print(f"coefficient arithmetic is done in {DefaultRing}")

    # display of various objects and operations (no self check [yet ?])

    # simple polynomials and basic arithmetic on them
    a = Polynomial([DefaultRing.elt(1), DefaultRing.elt(1)], DefaultRing)
    print(f"a = X + 1 = {a}")
    print(f"a * a = {a * a}")
    print(f"a - a = 0 = {a - a}")
    print(f"(a * a) % a = {(a * a) % a}")
    print(f"X^256 = {Polynomial.monomial(256, DefaultRing)}")
    print(f"mod = {mod}")
    print("f X^257 % mod = {Polynomial.monomial(257, DefaultRing) % mod}")

    # Random polynomial and basic NTT testing
    polys = [Polynomial.random(3, DefaultRing) for _ in range(2)]
    print(f"Random polys[:] = {list(str(p) for p in polys)}")
    ntts = [ntt_transform(poly) for poly in polys]
    print(f"ntt_transform(polys[0) = {ntts[0].coeffs}")
    print(f"fast_ntt_transform(polys[0) = {fast_ntt_transform(polys[0]).coeffs}")

    print(f"polys[0] + polys[1] = {polys[0] + polys[1]} (direct addition)" )
    print(f"polys[0] + polys[1] = {ntt_inv_transform(ntts[0] + ntts[1])} (ntt -> addition -> inv ntt)" )
    print(f"polys[0] * polys[1] % mod = {(polys[0] * polys[1]) % mod} (direct modulo multiplication)" )
    print(f"polys[0] * polys[1] % mod = {ntt_inv_transform(ntts[0] * ntts[1])} (ntt -> multiplication -> inv ntt)" )

    # fast NTT testing
    fast_ntts = [fast_ntt_transform(poly) for poly in polys]
    print(f"polys[0] * polys[1] % mod = {ntt_inv_transform(fast_ntts[0] * fast_ntts[1])} (fast ntt -> multiplication -> inv ntt)" )
    print(f"polys[0] * polys[1] % mod = {fast_inv_ntt_transform(fast_ntts[0] * fast_ntts[1])} (fast ntt -> multiplication -> fast inv ntt)" )