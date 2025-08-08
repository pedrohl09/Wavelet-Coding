import numpy as np
import math
import random
import time
from contextlib import ExitStack
import matplotlib.pyplot as plt  

class WaveletSystemSimulator:
    """
    Encapsula a simulação de um sistema de comunicação digital que usa
    códigos wavelet, modulação APK e um canal com desvanecimento Rayleigh.
    """

    def __init__(self, n_bits=100000, m=2, g=4):
        """
        Inicializa o simulador com os parâmetros do sistema.
        """
        # --- Parâmetros da Simulação ---
        self.n = n_bits
        self.m = m
        self.g = g
        self.mg = m * g
        self.M = 9

        # --- Matriz de Coeficientes Wavelet ---
        self.a = np.array([
            [1.0, 1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0],
            [1.0, 1.0, 1.0, -1.0, -1.0, -1.0, 1.0, -1.0]
        ])

        # --- Constelação APK (MPSK) ---
        self._create_mpsk_constellation()

        # --- Estado da Simulação ---
        simulation_length = self.n + self.mg - self.m
        self.x = np.zeros(self.n)
        self.y = np.zeros((simulation_length, self.m))
        self.simbr = np.zeros(simulation_length)
        
        self.w = np.zeros(self.m)
        self.z = np.zeros(self.mg)
        self.con = np.zeros(self.mg, dtype=int)
        self.idx = np.zeros(self.mg, dtype=int)
        
        # Variáveis de iteração
        self.iter_val = 0
        self.No = 0.0
        self.alfax = 0.0
        self.alfay = 0.0

    def _create_mpsk_constellation(self):
        PI = math.pi
        self.s = np.zeros(self.M); self.f = np.zeros(self.M); self.q = np.zeros(self.M); self.pb = np.zeros(self.M)
        self.s[0]=0; self.s[1]=2; self.s[2]=-2; self.s[3]=4; self.s[4]=-4
        self.f[0]=math.cos(0.0); self.f[1]=math.cos(0.316666667*PI); self.f[2]=math.cos(0.316666667*PI)
        self.f[3]=math.cos(0.55*PI); self.f[4]=math.cos(0.55*PI)
        self.q[0]=math.sin(0.0); self.q[1]=math.sin(0.316666667*PI); self.q[2]=-math.sin(0.316666667*PI)
        self.q[3]=math.sin(0.55*PI); self.q[4]=-math.sin(0.55*PI)
        self.s[5]=6; self.s[6]=-6; self.s[7]=8; self.s[8]=-8
        self.f[5]=math.cos(0.683333333*PI); self.f[6]=math.cos(0.683333333*PI); self.f[7]=math.cos(0.783333333*PI); self.f[8]=math.cos(0.783333333*PI)
        self.q[5]=math.sin(0.683333333*PI); self.q[6]=-math.sin(0.683333333*PI); self.q[7]=math.sin(0.783333333*PI); self.q[8]=-math.sin(0.783333333*PI)
        self.pb[0]=0.2734375; self.pb[1]=self.pb[2]=0.21875; self.pb[3]=self.pb[4]=0.109375; self.pb[5]=self.pb[6]=0.03125; self.pb[7]=self.pb[8]=0.00390625

    # ==============================================================================
    # MÉTODOS DE SIMULAÇÃO (Lógica Principal)
    # ==============================================================================

    def _fonte(self):
        return -1.0 if random.random() <= 0.5 else 1.0

    def _code(self):
        self.y[self.iter_val, :] = 0.0
        if self.iter_val < self.mg:
            for j in range(self.m):
                for i in range(self.iter_val // self.m + 1):
                    self.y[self.iter_val, j] += self.x[(i * self.m) + j] * self.a[j, self.iter_val - (self.m * i)]
        elif self.mg <= self.iter_val < self.n:
            res = self.iter_val % self.m
            for j in range(self.m):
                for i in range(self.g):
                    self.y[self.iter_val, j] += self.x[(self.iter_val - res - (self.m * i)) + j] * self.a[j, (self.m * i) + res]
        else:
            for j in range(self.m):
                for i in range(self.g - ((self.iter_val - self.n) // self.m) - 1):
                    self.y[self.iter_val, j] += self.x[(self.n - self.m) - (self.m * i) + j] * self.a[j, (self.iter_val - self.n) + (self.m * (i + 1))]

    def _modula(self, yenv):
        Amp=1.0; Ax,Ay=0.0,0.0
        for j in range(self.M):
            if yenv == self.s[j]:
                Ax=Amp*self.f[j]; Ay=Amp*self.q[j]; break
        return Ax,Ay

    def _soma_ruido(self, Sx, Sy):
        self.alfax=self._gauss(0,0.5); self.alfay=self._gauss(0,0.5)
        ruido_x=self._gauss(0,0.5*self.No); ruido_y=self._gauss(0,0.5*self.No)
        codx=Sx*self.alfax-Sy*self.alfay+ruido_x; cody=Sy*self.alfax+Sx*self.alfay+ruido_y
        return codx,cody

    def _demod(self, codx, cody):
        menor=float('inf')
        for j in range(self.M):
            est_x=self.f[j]*self.alfax-self.q[j]*self.alfay; est_y=self.q[j]*self.alfax+self.f[j]*self.alfay
            dist_sq=(codx-est_x)**2+(cody-est_y)**2; dist=(dist_sq/self.No)-math.log(self.pb[j])
            if dist<menor:
                self.simbr[self.iter_val]=self.s[j]; menor=dist
    
    def _decod(self):
        set_flag=0; yf=self.simbr[self.iter_val]
        for i in range(self.mg):
            res1=i%self.m; self.z[i]+=yf*self.a[res1,self.con[i]]*self.idx[i]; self.con[i]+=1
        o=0
        for i in range(self.mg):
            if self.con[i]==self.mg:
                self.w[o]=self.z[i]; set_flag=1; self.con[i]=0; self.z[i]=0; o+=1
        return set_flag

    def _gauss(self, mean, var):
        while True:
            v1=2.0*random.random()-1.0; v2=2.0*random.random()-1.0
            r=v1**2+v2**2
            if 0<r<1.0: break
        return mean+math.sqrt(var)*v1*math.sqrt(-2.0*math.log(r)/r)

    def run_simulation(self):
        random.seed(int(time.time() / 2))
        errdemod = np.zeros((81, 10))
        
        # --- 2. LISTAS PARA ARMAZENAR RESULTADOS PARA O GRÁFICO ---
        snr_values = []
        ber_values = []

        with ExitStack() as stack:
            try:
                fp_wav = stack.enter_context(open("wav_py_oop.txt", "w"))
                fp_hard = stack.enter_context(open("Hard2x8Ray_py_oop.txt", "w"))
            except IOError as e:
                print(f"Erro ao abrir arquivo: {e}")
                return

            print("Iniciando simulação... (SNR vs BER)")
            for snr in range(0, 31, 5):
                self.No = 10**(-0.1 * snr)
                vezes = 10 if snr <= 6 else 100
                erro_total = 0.0
                errdemsimbasimb = np.zeros((9, 9))

                for _ in range(vezes):
                    for i in range(self.n): self.x[i] = self._fonte()
                    for self.iter_val in range(self.n + self.mg - self.m): self._code()
                    for self.iter_val in range(self.n + self.mg - self.m):
                        yenv=np.sum(self.y[self.iter_val,:]); Ax,Ay=self._modula(yenv)
                        codx,cody=self._soma_ruido(Ax,Ay); self._demod(codx,cody)
                        simbaux=int(abs(self.simbr[self.iter_val])-1.0) if self.simbr[self.iter_val]<0 else int(self.simbr[self.iter_val])
                        yenvaux=int(abs(yenv)-1.0) if yenv<0 else int(yenv)
                        errdemsimbasimb[yenvaux,simbaux]+=1
                
                    pr=0; l=0; k=0; self.idx.fill(0); self.con.fill(0); self.z.fill(0)
                    for self.iter_val in range(self.n + self.mg - self.m):
                        if self.iter_val<=self.mg-self.m:
                            if self.iter_val==self.m*l:
                                l+=1
                                for _ in range(self.m): self.idx[k]=1; self.con[k]=0; k+=1
                        if self._decod():
                            for i in range(self.m):
                                if self.w[i]>=1.0: self.w[i]=1.0
                                elif self.w[i]<=-1.0: self.w[i]=-1.0
                                else: self.w[i]=self._fonte()
                                if 2*self.mg<=self.iter_val<self.n:
                                    if self.w[i]!=self.x[pr]: erro_total+=1
                                pr+=1

                total_bits_validos = vezes * (self.n - (2.0 * self.mg))
                ber = erro_total / total_bits_validos if total_bits_validos > 0 else 0
                
                # --- 3. ARMAZENANDO OS VALORES DE SNR E BER ---
                snr_values.append(snr)
                ber_values.append(ber)
                
                output_line = f"{float(snr):.4e} {ber:.4e}"
                fp_hard.write(output_line + "\n")
                print(output_line)
                fp_hard.flush()
            
            for i in range(81):
                fp_wav.write(" ".join(f"{val:.4e}" for val in errdemod[i, :8]) + "\n")
            
        print("\nSimulação concluída. Arquivos de saída e gráfico gerados.")

        # Plotagem dos resultados
        plt.figure(figsize=(10, 7))
        plt.semilogy(snr_values, ber_values, marker='o', linestyle='-', color='b', label='BER Simulada')
        plt.title('Desempenho do Sistema (BER vs. SNR)')
        plt.xlabel('SNR (dB)')
        plt.ylabel('Taxa de Erro de Bit (BER)')
        plt.grid(True, which="both", linestyle='--')
        plt.legend()
        plt.ylim(bottom=10e-7) # Limite inferior para melhor visualização
        plt.savefig("grafico_ber_vs_snr.png")
        print("\nGráfico salvo como 'grafico_ber_vs_snr.png'")

if __name__ == "__main__":
    simulator = WaveletSystemSimulator()
    simulator.run_simulation()