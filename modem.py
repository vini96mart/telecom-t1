import numpy as np


class Modem:
    def __init__(self, fs, bufsz, ans=False):
        self.bit_rate = 300
        self.bits = []
        self.fs = fs  # taxa de amostragem
        self.bufsz = bufsz  # quantidade de amostras que devem ser moduladas por vez

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
        if len(bits) == 0:
            bits.append(1)
        self.bits.extend(bits)


    def get_samples(self):
        x = np.arange(0, len(self.bits)/300, 1/self.fs)
        y = []
        t = 0.0
        phi = 0

        for i in self.bits:
            w = (self.rx_omega0 if i == 0 else self.rx_omega1)
            phi -= w*t
            for j in range(self.fs/self.bit_rate):
                x.append(t)
                y.append(np.sin(t+phi))
                t += 1/self.fs
        self.bits = []
        return np.array(y)

    # Demodulação

    def put_samples(self, data):
        for n in range (1, len(self.s)):
            self.v0r[n] = self.s[n] - self.r**self.L * np.cos(self.rx_omega0*self.L*self.T)*self.s[n-self.L] + self.r*np.cos(self.rx_omega0*self.T)*self.v0r[n-1] - self.r*np.sin(self.rx_omega0*self.T) * self.v0i[n-1]
            self.v0i[n] = -self.r**self.L * np.sin(self.rx_omega0*self.L*self.T)*self.s[n-self.L] + self.r*np.cos(self.rx_omega0*self.T)*self.v0i[n-1] - self.r*np.sin(self.rx_omega0*self.T) * self.v0r[n-1]

            self.v1r[n] = self.s[n] - self.r**self.L * np.cos(self.rx_omega1*self.L*self.T)*self.s[n-self.L] + self.r*np.cos(self.rx_omega1*self.T)*self.v1r[n-1] - self.r*np.sin(self.rx_omega1*self.T) * self.v1i[n-1]
            self.v1i[n] = -self.r**self.L * np.sin(self.rx_omega1*self.L*self.T)*self.s[n-self.L] + self.r*np.cos(self.rx_omega1*self.T)*self.v0i[n-1] - self.r*np.sin(self.rx_omega1*self.T) * self.v1r[n-1]
        rho = self.v1r**2+self.v1i**2+self.v0r**2+self.v0i**2
        pass

    def get_bits(self):
        return []
