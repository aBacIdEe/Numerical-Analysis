import math # for making up functions to integrate

class LegendrePolynomial:

    def __init__(self, coeff):
        self.coefficients = coeff
        self.roots = []

    def get_roots(self, prev_roots):
        if self.roots: return self.roots
        bounds = [-1] + prev_roots + [1]
        for i in range(len(bounds) - 1):
            self.roots.append(self.get_root((bounds[i], bounds[i + 1])))
        return self.roots

    def get_root(self, bound, last_guess=-2):
        # IVT
        if bound[0] == bound[1]: return bound[0]
        guess = (bound[0] + bound[1]) / 2
        if evaluate(self.coefficients, bound[0]) > evaluate(self.coefficients, bound[1]):
            higher, lower = bound[0], bound[1]
        else:
            higher, lower = bound[1], bound[0]

        value = evaluate(self.coefficients, guess)
        if guess == last_guess: return guess
        if (value) > 0:         return self.get_root((lower, guess), guess)
        elif (value) < 0:       return self.get_root((guess, higher), guess)
        else:                   return guess

def evaluate(coeffs, x):
    sum = 0
    for i, a in enumerate(coeffs):
        sum += a * (x ** (len(coeffs) - i - 1))
    return sum

def scalar_multiply(p, x):
    return [p[i] * x for i in range(len(p))]

def add_poly(p1, p2):
    if len(p1) >= len(p2):
        return [sum(n) for n in zip([0] * (len(p1) - len(p2)) + p2, p1)]
    else:
        return [sum(n) for n in zip([0] * (len(p2) - len(p1)) + p1, p2)]

legendrePolys = [LegendrePolynomial([1]), LegendrePolynomial([1, 0])]

def legendre_iterate():
    n = len(legendrePolys)
    temp1 = scalar_multiply(legendrePolys[-1].coefficients + [0], 2 * n - 1)
    temp2 = scalar_multiply(legendrePolys[-2].coefficients, -(n - 1))
    coeffs = scalar_multiply(add_poly(temp1, temp2), 1 / n)
    legendrePolys.append(LegendrePolynomial(coeffs))

def diff_poly(coeffs):
    if len(coeffs) <= 1: return [0]
    return [coeffs[i] * (len(coeffs) - i - 1) for i in range(len(coeffs))][:-1]

def get_roots(n):
    if n == 0: return []
    while len(legendrePolys) - 1 < n:
        legendre_iterate()
    return legendrePolys[n].get_roots(get_roots(n - 1))

def get_weights(n):
    ret, roots = [], get_roots(n)
    for root in roots:
        temp = evaluate(diff_poly(legendrePolys[n].coefficients), root) ** 2
        ret.append(2 / ((1 - root ** 2) * temp))
    return ret

def gaussQ(f, a, b, w, x):
    sum = 0
    for i in range(len(w)):
        sum += ((b - a) / 2) * f((b - a) / 2 * x[i] + (a + b) / 2) * w[i]
    return sum

# the function you want to integrate
def func(x):
    return math.log(x)

print(gaussQ(func, 0, math.pi, get_weights(40), get_roots(40)))
print(get_weights(10))