from fractions import Fraction

fact = [1]  # 階乗を保持


def fac(n):  # 階乗を計算
    while len(fact) <= n:
        fact.append(fact[-1] * len(fact))
    return fact[n]


def com(n, k):  # 二項係数計算パート
    return fac(n) // (fac(k) * fac(n - k))


def add(a, b):  # 多項式の足し算 O(|a| + |b|)
    c = [Fraction(0) for _ in range(max(len(a), len(b)))]
    for i in range(len(a)):
        c[i] += a[i]
    for i in range(len(b)):
        c[i] += b[i]
    return c


def sub(a, b):  # 多項式の引き算 O(|a| + |b|)
    c = [Fraction(0) for _ in range(max(len(a), len(b)))]
    for i in range(len(a)):
        c[i] += a[i]
    for i in range(len(b)):
        c[i] -= b[i]
    return c


def mul(a, b):  # 多項式の掛け算 O(|a| * |b|)
    c = [Fraction(0) for _ in range(len(a) + len(b) - 1)]
    for i in range(len(a)):
        for j in range(len(b)):
            c[i + j] += a[i] * b[j]
    return c


N, M = map(int, input().split())  # 入力の受け取り \int_0^1 x^{N+2M}/artanh^N(x) dx
M2 = M << 1
N2M = N + M2
N2M2 = N2M << 1  # 次数の管理 (この程度になる)

DI = [[Fraction(0) for i in range(N2M)] for j in range(N)]  # 部分積分パート
for k in range(M2 + 1):  # 初期値の設定
    DI[0][k] = Fraction(com(M2, k))
    if k & 1:
        DI[0][k] *= -1
for _ in range(N - 1):  # O(N^2 * (N + M))
    DD = [[Fraction(0) for i in range(N2M)] for j in range(N)]
    for k in range(N):
        for m in range(N2M):
            if DI[k][m] == 0:
                continue
            kei = DI[k][m] / (k + 1)
            DD[0][m + 1] += kei
            DD[k + 1][m + 1] -= kei
            DD[k + 1][m] += kei
    DI = DD

S_L = [[1] for i in range(N2M)]  # ガンマ関数の正規化のための左右からの累積積
for i in range(1, N2M):  # O((N + M)^2)
    S_L[i] = mul(S_L[i - 1], [i, 1])
S_R = [[1] for i in range(N2M)]
for i in reversed(range(0, N2M - 1)):  # O((N + M)^2)
    S_R[i] = mul(S_R[i + 1], [N2M - i, -1])
s = [0]  # \int_0^1 s * Gamma(1+x)Gamma(2-x) dx
for i2 in range(N2M):  # O((N + M)^3)
    ad = [DI[i1][i2] for i1 in range(N)]
    s = add(s, mul(ad, mul(S_L[i2], S_R[i2])))

II = [[] for _ in range(N2M2)]  # \int_0^{pi/2} x^n * cos2(2k+1)xdx の k を並べたもの
JJ = [[] for _ in range(N2M2)]  # \int_0^{pi/2} x^n * sin2(2k+1)xdx の k を並べたもの
II[0] = [Fraction(0), Fraction(0)]
JJ[0] = [Fraction(0), Fraction(1)]
for i in range(1, N2M2):  # O((N + M)^2)
    II[i] = [Fraction(0) for _ in range(i + 2)]
    JJ[i] = [Fraction(0) for _ in range(i + 2)]
    for j in range(i + 1):
        II[i][j + 1] -= JJ[i - 1][j]
        JJ[i][j + 1] += II[i - 1][j]
    JJ[i][1] += Fraction(1, 2 * fac(i))
B = [[] for _ in range(N2M2)]
for i in range(1, N2M2):  # O((N + M)^2)
    B[i] = [Fraction(0) for _ in range((i >> 1) + 1)]
    MUL = (i + 1) * 4 * fac(i)
    for j in range(2, i + 2, 2):
        B[i][(j - 2) >> 1] = - MUL * (1 - Fraction(1, 1 << (j + 1))) * II[i][j]
ans = []
for i in range(len(s)):
    ad = sub(B[i + 1], B[i])
    for idx in range(len(ad)):
        ad[idx] *= s[i]
    ans = add(ans, ad)
last = Fraction(1 << N, fac(N2M + 1))
for i in range(len(ans)):
    ans[i] *= last
# output_1
out1 = ""
for i in reversed(range(len(ans))):
    if ans[i] == 0:
        continue
    if ans[i].numerator > 0:
        out1 += '+'
    if ans[i].denominator == 1:
        out1 += str(ans[i].numerator)
    else:
        bu = ans[i].numerator
        if bu < 0:
            out1 += '-'
            bu = -bu
        out1 += '\\frac{' + str(bu) + '}{' + str(ans[i].denominator) + '}'
    out1 += '\\frac{\\zeta(' + str(2 * i + 3) + \
        ')}{\\pi^{' + str(2 * i + 2) + '}}'
print(out1)
# output_2
out2 = ""
for i in reversed(range(len(ans))):
    if ans[i] == 0:
        continue
    if ans[i].numerator > 0:
        out2 += '+'
    if ans[i].denominator == 1:
        out2 += str(ans[i].numerator) + '*'
    else:
        bu = ans[i].numerator
        if bu < 0:
            out2 += '-'
            bu = -bu
        out2 += str(bu) + '/' + str(ans[i].denominator) + '*'
    out2 += 'zeta(' + str(2 * i + 3) + ')/pi^(' + str(2 * i + 2) + ')'
print(out2)
