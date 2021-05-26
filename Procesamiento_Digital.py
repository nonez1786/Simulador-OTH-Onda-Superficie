# -*- coding: utf-8 -*-
"""
Created on Thu May 13 23:31:45 2021

@author: Zenon
"""
#import commpy
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftfreq
from numpy import r_
import Generador_Senal




# Transmisor  
Fc = 0 # [Hz]
AB = 20e3 # [Hz]
T_pulso = 800e-6 # [s]
N = 2 # Num de integraciones
num_samples =1 #int(T_pulso*fs_banda_base)
f_start = 0 # [Hz]
f_stop = AB # [Hz]
PRF = 0 # [Hz]




#=============================================================================

def rrcosfilter(N, alpha, Ts, Fs):
    """
    Generates a root raised cosine (RRC) filter (FIR) impulse response.
    Parameters
    ----------
    N : int
        Length of the filter in samples.
    alpha : float
        Roll off factor (Valid values are [0, 1]).
    Ts : float
        Symbol period in seconds.
    Fs : float
        Sampling Rate in Hz.
    Returns
    ---------
    time_idx : 1-D ndarray of floats
        Array containing the time indices, in seconds, for
        the impulse response.
    h_rrc : 1-D ndarray of floats
        Impulse response of the root raised cosine filter.
    """

    T_delta = 1/float(Fs)
    time_idx = ((np.arange(N)-N/2))*T_delta
    sample_num = np.arange(N)
    h_rrc = np.zeros(N, dtype=float)

    for x in sample_num:
        t = (x-N/2)*T_delta
        if t == 0.0:
            h_rrc[x] = 1.0 - alpha + (4*alpha/np.pi)
        elif alpha != 0 and t == Ts/(4*alpha):
            h_rrc[x] = (alpha/np.sqrt(2))*(((1+2/np.pi)* \
                    (np.sin(np.pi/(4*alpha)))) + ((1-2/np.pi)*(np.cos(np.pi/(4*alpha)))))
        elif alpha != 0 and t == -Ts/(4*alpha):
            h_rrc[x] = (alpha/np.sqrt(2))*(((1+2/np.pi)* \
                    (np.sin(np.pi/(4*alpha)))) + ((1-2/np.pi)*(np.cos(np.pi/(4*alpha)))))
        else:
            h_rrc[x] = (np.sin(np.pi*t*(1-alpha)/Ts) +  \
                    4*alpha*(t/Ts)*np.cos(np.pi*t*(1+alpha)/Ts))/ \
                    (np.pi*t*(1-(4*alpha*t/Ts)*(4*alpha*t/Ts))/Ts)

    return time_idx, h_rrc


#=============================================================================

def correlacion(senal_1,senal_2): 
    """
    Realiza la funcion correlacion entre las senal_1 y senal_2
    
    array correlacion(array senal_1,array senal_2)
    Entradas: 
             senal_1: vector de longitusd N
             senal_2: vector de longitud M
    Salidas: 
             salida: vector que contiene la correlacion entre las señales de entradas. 
             Tiene una longitud N + M -1
    """
    salida = np.convolve(senal_1,np.flip(senal_2.conjugate()))
    return salida


#=============================================================================
def FFT(y,n,dt):
    """
    Realiza la FFT de la señal y
    
    array(complex) FFT(array complex y,int n, float dt)
    Entradas: 
             y: vector de longitud N
             n: Numero de muestras del vector y
             dt: intevalo de muestreo de la señal y
    Salidas: 
             Y : resultado de FFT de y
             frec: Frecuencias asociadas a FFT(y) 
    """
    Y = fft(y) / n # Normalizada
    frec = fftfreq(n, dt) # Recuperamos las frecuencias
    return Y,frec
  
    
def graficar_t(t,y,imag_real):
    plt.plot(t*1e3,y,label = imag_real)
    plt.xlabel('Tiempo (ms)'),plt.ylabel('Amplitud'),plt.grid()
    plt.legend()
    
def graficar_f(f,Y,Fc,AB):
    plt.vlines(f/1000, 0, np.abs(Y),'r')
    #plt.ylim(-0.1, 1),
    plt.xlim(-4*(Fc+AB)/1000,4*(Fc+AB)/1000)
    plt.xlabel('Frecuencia (kHz)'),plt.ylabel('Abs(Y)'),plt.grid()
    plt.vlines((Fc+AB/2)/1000, 0, 0.5,'g') 
    plt.vlines((Fc-AB/2)/1000, 0, 0.5,'g')     
    



#=============================================================================
# Funcion Principal
def main():
    print('\n=============================================================')
    print('                        INICIO                                 ')
    print('=============================================================\n')
    
    
    [senal1,t] = Generador_Senal.Generacion_Tx("LFM","None",N,AB,PRF,T_pulso,Fc)
    #[senal1,t] = Generador_Senal.Generacion_Tx("Pulse_Cod","Complem_1",N,AB,PRF,T_pulso,Fc)
    
        
    ts = t[1]-t[0]
    plt.figure(1,figsize = (10,5))
    plt.subplot(211)   
    graficar_t(t,senal1.real,"Señal BPSK")
       
    plt.subplot(212)  
    [Y,frq] = FFT(senal1,np.size(t),ts)
    graficar_f(frq,Y,Fc,AB)
    
    
    
    """
    
    filtro = rrcosfilter(1*np.size(senal1), 0.35,1/AB,1/ts)
    plt.figure(2,figsize = (10,5))
    plt.subplot(211)  
    graficar_t(filtro[0],filtro[1], "Filtro RRC")
    
    plt.subplot(212)  
    [Y,frq] = FFT(filtro[1],np.size(filtro[0]),ts)
    graficar_f(frq, Y*5, Fc, AB)
    """    
    
    """
    plt.figure(3,figsize = (10,5))
    y = np.convolve(filtro[1],senal2,mode='same')
    t_filtro = np.arange(0,np.size(y)*ts,ts)
    
    graficar_t(t_filtro, y, "Pulso_RRC")
    graficar_t(t,senal2.real*10,"Pulso_Rec")
    """
    
    """
    plt.figure(3,figsize = (10,5))
    [senal_tx,tx] = Generador_Senal.Generacion_Tx("Pulse_Cod","Complem_1",1,AB,0,T_pulso,Fc)
    Salida_Filtro_Adaptax = correlacion(senal1, senal_tx)
    Salidatx = np.linspace(0, (np.size(Salida_Filtro_Adaptax) - 1) * (ts), np.size(Salida_Filtro_Adaptax))
    graficar_t(Salidatx, Salida_Filtro_Adaptax, "Auto_Corre_Cod_1")
    
    """    
    
    
    
# INICIO    
if __name__ == '__main__':
    main()



