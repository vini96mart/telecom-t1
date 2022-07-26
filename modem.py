import numpy as np
import scipy.signal


class Modem:
    def __init__(self, fs, bufsz, ans=False):
        self.bit_rate = 300
        self.bits = []
        self.s = []
        self.fi = 0
        self.sbuffer = []
        self.lst = []

        self.v0r = np.zeros(bufsz)
        self.v0i = np.zeros(bufsz)
        self.v1r = np.zeros(bufsz)
        self.v1i = np.zeros(bufsz)
        self.y = []
        self.v = []
        self.fs = fs  # taxa de amostragem
        self.bufsz = bufsz  # quantidade de amostas que devem ser moduladas por vez

        # frequências de modulação (upload)
        self.tx_omega0 = 2*np.pi*(1080 + 100)
        self.tx_omega1 = 2*np.pi*(1080 - 100)
        # frequências de demodulação (download)
        self.rx_omega0 = 2*np.pi*(1750 + 100)
        self.rx_omega1 = 2*np.pi*(1750 - 100)
        # se o modem estiver atendendo uma ligação
        if ans:
            # inverte as frequências
            self.tx_omega0, self.rx_omega0 = self.rx_omega0, self.tx_omega0
            self.tx_omega1, self.rx_omega1 = self.rx_omega1, self.tx_omega1

    # Modulação

    def put_bits(self, bits):
        self.bits.extend(bits)

    def get_samples(self):
        y = []

        if len(self.bits) == 0:
            self.bits.append(1)

        if self.bits[0] != 0:
            w = self.tx_omega1
        else:
            w = self.tx_omega0

        for j in range(self.bufsz):
            y.append(np.sin(self.fi))
            self.fi += w/self.fs

        self.bits = self.bits[1:]
        return np.array(y)
    # Demodulação

    def put_samples(self, data):
        self.s.extend(data)

    def get_bits(self):
        fs = self.fs
        s = []
        aux = [0] * self.bufsz
        s += self.sbuffer
        if len(s) == 0:
            s += aux
        s.extend(self.s)
        self.sbuffer = s[len(s)-self.bufsz:]

        T = 1/fs
        L = fs//300
        omega0, omega1 = self.rx_omega0, self.rx_omega1
        r = 0.99
        rr = 0.9999

        rL = r**L

        c0LT, s0LT, c0T, s0T = np.cos(
            omega0*L*T), np.sin(omega0*L*T), np.cos(omega0*T), np.sin(omega0*T)
        c1LT, s1LT, c1T, s1T = np.cos(
            omega1*L*T), np.sin(omega1*L*T), np.cos(omega1*T), np.sin(omega1*T)

        size = len(self.v0r)

        self.v0i = np.append(self.v0i, np.zeros(len(self.s)))
        self.v0r = np.append(self.v0r, np.zeros(len(self.s)))
        self.v1i = np.append(self.v1i, np.zeros(len(self.s)))
        self.v1r = np.append(self.v1r, np.zeros(len(self.s)))

        v0r, v0i, v1r, v1i = self.v0r, self.v0i, self.v1r, self.v1i

        for n in range(size, len(s)):
            v0r[n] = s[n] - r**L * c0LT*s[n-L] + r * \
                np.cos(omega0*T)*v0r[n-1] - r*np.sin(omega0*T)*v0i[n-1]
            v0i[n] = -rL * s0LT * s[n-L] + r*c0T*v0i[n-1] + r * s0T * v0r[n-1]

            v1r[n] = s[n] - rL * c1LT * s[n-L] + \
                r * c1T*v1r[n-1] - r * s1T * v1i[n-1]
            v1i[n] = -rL * s1LT * s[n-L] + r * \
                c1T * v1i[n-1] + r * s1T * v1r[n-1]

        v0 = v0i**2 + v0r**2
        v1 = v1i**2 + v1r**2

        c = abs(v1 - v0)
        v = np.zeros(len(c))
        y = np.zeros(len(c))

        if len(self.v) > 0:
            v[0], v[1] = self.v[0], self.v[1]

        if len(self.y) > 0:
            y[0], y[1] = self.y[0], self.y[1]

        for n in range(2, len(c)):
            v[n] = (1-rr) * c[n] + 2 * rr * \
                np.cos(2*np.pi*300/fs) * v[n-1] - rr**2 * v[n-2]
            y[n] = v[n] - v[n-2]

        self.v = v[len(v)-2:]
        self.y = y[len(y)-2:]

        delta = v1 - v0

        point = (y[1:] > 0) & (y[:-1] < 0)
        bits = []

        for i in range(len(point)):
            if point[i] != 0 and not self.lst.count(i):
                if delta[i] > 0:
                    bits.append(1)
                else:
                    bits.append(0)
                for j in range(-self.bufsz//4, self.bufsz//4+1):
                    self.lst.append(i+j)

        for i in range(len(self.lst)):
            self.lst[i] = self.lst[i] - self.bufsz

        if(len(self.lst) > 0):
            while(len(self.lst) > 0 & self.lst[0] < 0):
                self.lst.pop(0)

        self.s = []
        self.v0i = v0i[len(v0i)-self.bufsz:]
        self.v0r = v0r[len(v0r)-self.bufsz:]
        self.v1i = v1i[len(v1i)-self.bufsz:]
        self.v1r = v1r[len(v1r)-self.bufsz:]

        return bits
