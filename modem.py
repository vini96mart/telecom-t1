import numpy as np
import scipy.signal

class Modem:
    def __init__(self, fs, bufsz, ans=False):
        self.bit_rate = 300
        self.bits = []
        self.s = []
        self.phi = 0
        self.fs = fs  # taxa de amostragem
        self.bufsz = bufsz  # quantidade de amostras que devem ser moduladas por vez

        self.v0r = np.zeros(len(self.s))
        self.v0i = np.zeros(len(self.s))
        self.v1r = np.zeros(len(self.s))
        self.v1i = np.zeros(len(self.s))
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
        

        w = (self.rx_omega0 if self.bits[0] == 0 else self.rx_omega1)
        for j in range(int(self.fs/self.bit_rate)):
            y.append(np.sin(self.phi))
            self.phi += w/self.fs
        self.bits = self.bits[1:]
        return np.array(y)

    # Demodulação

    def put_samples(self, data):
        self.s.extend(data)

    def get_bits(self):
        fs = self.fs
        s = self.s
        T = 1/fs
        L = fs//300
        omega0, omega1 = self.rx_omega0, self.rx_omega1
        r = 0.99
        rr = 0.9999
        
        v0r, v0i, v1r, v1i = self.v0r, self.v0i, self.v1r, self.v1i
        rL = r**L
        
        c0LT, s0LT, c0T, s0T = np.cos(omega0*L*T), np.sin(omega0*L*T), np.cos(omega0*T), np.sin(omega0*T)
        c1LT, s1LT, c1T, s1T = np.cos(omega1*L*T), np.sin(omega1*L*T), np.cos(omega1*T), np.sin(omega1*T)

        for n in range (L, len(self.s)):
            v0r = s[n] - rL*c0LT*s[n-L] + r*c0T*v0r - r*s0T*v0i
            v0i = -rL*s0LT*s[n-L] + r*c0T*v0i - r*s0T*v0r

            v1r = s[n] - rL *c1LT*s[n-L] + r*c1T*v1r - r*s1T* v1i
            v1i = -rL *s1LT*s[n-L] + r*c1T*v0i - r*s1T*v1r

        v0 = v0r**2 + v0i**2
        v1 = v1r**2 + v1i**2

        rho = v0 + v1

        c = abs(v1 - v0)
        v = np.zeros(len(c))
        y = np.zeros(len(c))

        for n in range (1, len(s)):
            v[n] = (1-rr)*c[n] + 2*rr*np.cos(2*np.pi*300/fs)*v[n-1] - rr**2*v[n-2]
            y[n] = v[n] - v[n-2]

        delta = rho

        filt = scipy.signal.firwin(40, 300, pass_zero='lowpass', fs = self.fs)
        ponto_amostra = 1500*((y[1:] > 0) & (y[:-1] < 0))
        bits = []

        for i in range (len(ponto_amostra)):
            if ponto_amostra[i] != 0:
                bits.append(1 if delta[i] > 0 else 0)
        
        self.s = []
        
        #print("bits = ", bits)
        #plt.show()
        return bits
