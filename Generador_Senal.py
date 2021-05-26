#==============================================================================
#Titulo: Generador_Senal
#Autor: Zenon Saavedra 
#==============================================================================
"""
En este se encuentran todas funciones encargadas de la generacion de las señales:
    Eco
    Ruido
    Clutter
    Interferencia
    Recibida
    
Se ingresan un conjunto de parametros,y como salida se tienen series Temporales

Created on Fri Mar 19 10:46:53 2021
@author: Zenon Saavedra
"""

import os
import sys
import numpy as np
from numpy import r_
import math
from math import pi
from numpy import sin
from numpy import cos



# Codigos
Cod_complem_1 = np.array([1,1,0,1,1,1,1,0,1,0,0,0,1,0,1,1])
Cod_complem_2 = np.array([1,1,0,1,1,1,1,0,0,1,1,1,0,1,0,0])
Cod_Barker = np.array([1,1,1,0,0,1,0])
Cod_Barker_7 = np.array([1,1,1,0,0,0,1,0,0,1,0])
Cod_Barker_13 = np.array([1,1,1,1,1,0,0,1,1,0,1,0,1])
Cod_Frank_1 = np.array([1,1,1,-1])
Cod_Frank_2 = np.array([1,1,1,1,-0.5+0.87j,-0.5-0.87j,1,-0.5-0.87j,-0.5+0.87j])
Cod_Frank_3 = np.array([1,1,1,1,1,1j,-1,-1j,1,-1,1,-1,1,-1j,-1,1j])


#------------------------------------------------------------------------------
def Generacion_Tx(Tipo_radar,Tipo_cod,N,AB,PRF,T,fc):
    """
        Genera la serie temporal Tx, de acuerdo al Tipo de Radar (CW o Pulsado)
            
        array complex Generacion_Tx( string Tipo_radar, string Tipo_cod, int N, float AB, float PRF, float T, float fc)
        Entrada:
            Tipo_radar: "Pulsado"  "CW"
            Tipo_cod:"Complementario"  "Barker"  "Frank "
            N: Num intergraciones
            AB: Ancho de Banda
            PRF: Frecuencia de Repeticion de Pulso
            T: Ancho del pulso
            fc: Frecuencia de portadora
        Salida: 
            senal: muestras de la Serie Temporal      
            t: base de tiempo, arrancando con un determinado offset
    """   
    #global senal
    senal = []
    t = []
    ofset = 0
    if(Tipo_radar == 'LFM'):
        fase = 0
        [aux, auxt] = Generacion_LFM(fc, T,  AB, PRF, fase)
        for i in range(N):
            senal = np.hstack((senal,aux))
            #t = np.hstack((t,auxt))
            t = np.append(np.array(t), auxt + ofset)
            ofset = np.array(t)[-1]
            
    elif(Tipo_radar == 'Pulse_Cod'):           
        [aux, auxt]= Generacion_Pulso(Tipo_cod, fc, AB, T, PRF)
        for i in range(N):
            senal = np.hstack((senal,aux))
            #t = np.hstack((t,auxt))
            t = np.append(np.array(t), auxt + ofset)
            ofset = np.array(t)[-1]
    else:
        print("\n Error: opcion No valida. Tipo de Radar\n") 
              
    return senal,t
#------------------------------------------------------------------------------

def Generacion_LFM(fc, T, AB, PRF, fase):
    """
        Genera una serie temporal, correspondiente a un señal LFM-UP
            
        array complex Generacion_LFM(float fc, float T, float AB, int fase)
        Entrada:
            fc: Frecuencia de portadora [Hz]
            T: Ancho del pulso [s]
            AB: Ancho de Banda [Hz]
            PRF :Frec de Preteticon de Pulso [Hz]
            fase: Fase Inicial [rad]
        Salida: 
            senal: muestras de la Serie Temporal LFM      
            t: base de tiempo, arrancando con un determinado offset
    """ 
    
    fmax = fc
    if (fc == 0):
        fmax = AB
    fs = 16*fmax
    num_muestras = int (T*fs)
        
    f_inicial = 0 
    f_final = AB
    times_s = np.linspace(0, T, num_muestras) # Chirp tiempo.
    k = (f_final - f_inicial) / T # Chirp rate.
    sweepFreqs_Hz = (fc + f_inicial + k * times_s) * times_s
    senal = np.exp(1j* (fase + 2 * np.pi * sweepFreqs_Hz))

    t = np.linspace(0, (num_muestras - 1) * (1/fs), num_muestras) # Intervalo de tiempo en segundos
    
    if (PRF != 0):
        PRP = 1/PRF
        ts = 1/fs
        t_recepcion = np.arange(t[-1],int(PRP/ts)*ts,ts)
        recepcion = np.zeros(np.size(t_recepcion))
        t = np.append(t,t_recepcion)
        senal = np.append(senal,recepcion)
    

    return senal,t
#------------------------------------------------------------------------------
 


def Barker(fc,AB,PRF):
    print('Barker 7')
    Cod_Barker = np.array([[1],[1],[1],[0],[0],[1],[0]]) # Araay (7,1) no es array 1D
    senal = senal_BPSK(Cod_Barker,fc,AB,PRF)
    return senal  

def Barker7(fc,AB,PRF):
    print('Barker 11')
    Cod_Barker_7 = np.array([[1],[1],[1],[0],[0],[0],[1],[0],[0],[1],[0]])
    senal = senal_BPSK(Cod_Barker_7,fc,AB,PRF)
    return senal  

def Barker13(fc,AB,PRF):
    print('Barker 13')
    Cod_Barker_13 = np.array([[1],[1],[1],[1],[1],[0],[0],[1],[1],[0],[1],[0],[1]])
    senal = senal_BPSK(Cod_Barker_13,fc,AB,PRF)
    return senal  
    
def Complem_1(fc,AB,PRF):
    print('Complementario_1')
    Cod_complem_1 = np.array([[1],[1],[0],[1],[1],[1],[1],[0],[1],[0],[0],[0],[1],[0],[1],[1]])
    senal = senal_BPSK(Cod_complem_1,fc,AB,PRF)
    return senal
    
def Complem_2(fc,AB,PRF):
    print('Complementario_2')                            
    Cod_complem_2 = np.array([[1],[1],[0],[1],[1],[1],[1],[0],[0],[1],[1],[1],[0],[1],[0],[0]])
    senal = senal_BPSK(Cod_complem_2,fc,AB,PRF)
    return senal

def Frank1(fc,AB,PRF):
    print('Frank 1')
    return senal

def Frank2(fc,AB,PRF):
    print('Frank 2')
    return senal

def Frank3(fc,AB,PRF):
    print('Frank 3')
    return senal



#------------------------------------------------------------------------------
    
def Generacion_Pulso(Tipo_cod,fc,AB,T,PRF):
    """
        Se genera señales BPSK codificadas, en funcion del codigo elegido.
        
        array complex Generacion_Pulso(string Tipo_Cod, float fc, float AB, float T, float PRF)
        Entrada:
            Tipo_Cod: Complemetario/Barker7,11,13/Frank 
            fc: Frecuencia de Portadora [Hz]
            AB: Ancho de Banda [Hz]
            T: Ancho de Pulso [s]
            PRF: Frecuencia de Repeticion de Pulso [Hz]
        Salida: 
            senal: muestras de la Serie Temporal PBSK Codificada       
            t: base de tiempo
    """

    senal = []    
    if(Tipo_cod == 'Complem_1'):
        [senal, t] = Complem_1(fc,AB,PRF)
    
    elif(Tipo_cod == 'Complem_2'):
        [senal, t] = Complem_2(fc,AB,PRF)
        
    elif(Tipo_cod == 'Barker_7'):           
        [senal, t] = Barker(fc,AB,PRF)
        
    elif(Tipo_cod == 'Barker_11'):           
        [senal, t] = Barker7(fc,AB,PRF)
        
    elif(Tipo_cod == 'Barker_13'):           
        [senal, t] = Barker13(fc,AB,PRF)
        
    else:
        print("\n Error: opcion No valida. Tipo de Codigo\n") 

    return senal,t
#------------------------------------------------------------------------------

def senal_BPSK(codigo,fc,AB,PRF):
    """
        Genera una serie temporal, correspondiente a un señal BPSK-Codificda
        
        array complex senal_BPSK(array int codigo, float fc, float AB)

        Entrada:
            codigo: Codigo modulador
            fc: Frecuencia de portadora [Hz]
            AB: Ancho de Banda [Hz]
            PRF: Frecuencia de repeicion de pulso [Hz]
        Salida: 
            senal: muestras de la Serie Temporal PBSK Codificada       
            tiempo: base de tiempo, arrancando con un determinado offset
    """
    Nbits = np.size(codigo,0) # num bits
    T = Nbits*(1/AB)
    fmax = fc
    if (fc == 0):
        fmax = AB        
        
    fs = 16*fmax
    N = int(T*fs) # Num de muestras del pulso
    Ns = int(fs/AB)  # Num de muestras por bit
    bits = codigo>0
    M = np.tile(bits*2-1,(1,Ns))
    t = np.linspace(0, (N - 1) * (1/fs), N) # t = r_[0.0:N]/fs
    bpsk_I = M.ravel()*cos(2*pi*fc*t)
    bpsk_C = M.ravel()*sin(2*pi*fc*t)
    senal = bpsk_I + (1j)*bpsk_C
    
    if (PRF != 0):
        PRP = 1/PRF
        ts = 1/fs
        t_recepcion = np.arange(t[-1],int(PRP/ts)*ts,ts)
        recepcion = np.zeros(np.size(t_recepcion))
        t = np.append(t,t_recepcion)
        senal = np.append(senal,recepcion)
    
    return senal,t
       
#------------------------------------------------------------------------------
def ruido():
    print("Ruido")
    return ruido
    
#------------------------------------------------------------------------------
def clutter():
    print("Clutter")
    return clutter

#------------------------------------------------------------------------------
def Eco():
    print("Eco")
    return Eco

def Senal_Recibida():
    print("Señal Recibida")
    return Senal_Rx
    



