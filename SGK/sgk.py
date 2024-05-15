import pandas as pd
from findpeaks import findpeaks
import numpy as np
import os  
import time
import matplotlib.pyplot as plt
import sys
#from astropy.stats import LombScargle
from astropy.timeseries import LombScargle
from PyAstronomy.pyTiming import pyPDM
import libwwz
import libwwz.wwz as wwz

from ipywidgets import FloatSlider
from ipywidgets import BoundedFloatText
from ipywidgets import interact

from SGK.gwwz import WWZ
from SGK.gwwzplotter import WWZPlotter

from scipy.optimize import curve_fit



def sumar(a,b):
    return a+b
################################################################################################################################################################################################################################################

def lista_peak_fp(pp, tt):
    fp = findpeaks(lookahead=1)

    results = fp.fit(pp)

    peak_x = results['df']['x']
    peak_y = results['df']['y']
    peak_values = results['df']['peak']

    all_peaks = [(tt[x], y, value) for x, y, value in zip(peak_x, peak_y, peak_values)]

    filtered_rows = [row for row in all_peaks if row[2] == True]
    
    sorted_rows = sorted(filtered_rows, key=lambda x: x[1], reverse=True)

    modified_data_list = [[1/row[0], row[0], row[1]] for row in sorted_rows]

    df = pd.DataFrame(modified_data_list, columns=['Periodos', 'Freq', 'Power'])

    print("DataFrame:")
    print(df)
    

    fp.plot()

    return df


def lista_peak(pp, tt):
    fp = findpeaks(lookahead=1)

    results = fp.fit(pp)

    peak_x = results['df']['x']
    peak_y = results['df']['y']
    peak_values = results['df']['peak']

    all_peaks = [(tt[x], y, value) for x, y, value in zip(peak_x, peak_y, peak_values)]

    filtered_rows = [row for row in all_peaks if row[2] == True]
    
    sorted_rows = sorted(filtered_rows, key=lambda x: x[1], reverse=True)

    modified_data_list = [[ row[0], row[1]] for row in sorted_rows]

    df = pd.DataFrame(modified_data_list, columns=['Axis X', 'Axis Y'])

    print("DataFrame:")
    print(df)
    

    fp.plot()

    return df



####################################################################################################################################################################################

###def saca_punto(ff,pp,N):
    sys.stdout = open("nul", "w")

    fp = findpeaks(lookahead=1)
    results = fp.fit(pp)
    
    peak_x = results['df']['x']
    peak_y = results['df']['y']
    peak_values = results['df']['peak']
    
    all_peaks = [(ff[x], y, value) for x, y, value in zip(peak_x, peak_y, peak_values)]
    filtered_rows = [row for row in all_peaks if row[2] == True]
    sorted_rows = sorted(filtered_rows, key=lambda x: x[1], reverse=True)
    sorted_rows = sorted_rows[:N]
    modified_data_list = [[row[0], row[1]] for row in sorted_rows]
    sys.stdout = sys.__stdout__
    return modified_data_list

###def saca_punto2(ff,pp):
    sys.stdout = open("nul", "w")

    fp = findpeaks(lookahead=1)
    results = fp.fit(pp)
    
    peak_x = results['df']['x']
    peak_y = results['df']['y']
    peak_values = results['df']['peak']
    
    all_peaks = [(ff[x], y, value) for x, y, value in zip(peak_x, peak_y, peak_values)]
    filtered_rows = [row for row in all_peaks if row[2] == True]
    sorted_rows = sorted(filtered_rows, key=lambda x: x[1], reverse=True)
    modified_data_list = [[row[0], row[1]] for row in sorted_rows]
    sys.stdout = sys.__stdout__
    return modified_data_list

def saca_punto(ff,pp,N):

    fp = findpeaks(lookahead=1)
    results = fp.fit(pp)
    
    peak_x = results['df']['x']
    peak_y = results['df']['y']
    peak_values = results['df']['peak']
    
    all_peaks = [(ff[x], y, value) for x, y, value in zip(peak_x, peak_y, peak_values)]
    filtered_rows = [row for row in all_peaks if row[2] == True]
    sorted_rows = sorted(filtered_rows, key=lambda x: x[1], reverse=True)
    sorted_rows = sorted_rows[:N]
    modified_data_list = [[row[0], row[1]] for row in sorted_rows]

    return modified_data_list

def saca_punto2(ff,pp):

    fp = findpeaks(lookahead=1)
    results = fp.fit(pp)
    
    peak_x = results['df']['x']
    peak_y = results['df']['y']
    peak_values = results['df']['peak']
    
    all_peaks = [(ff[x], y, value) for x, y, value in zip(peak_x, peak_y, peak_values)]
    filtered_rows = [row for row in all_peaks if row[2] == True]
    sorted_rows = sorted(filtered_rows, key=lambda x: x[1], reverse=True)
    modified_data_list = [[row[0], row[1]] for row in sorted_rows]
    return modified_data_list


def alaWWZ(WWZ_simple_linear,N):
    l1 = []
    l2 = []
    l3 = []
    l4 = []
    l5 = []
    for i in range(len(WWZ_simple_linear[0])):
        ff = WWZ_simple_linear[1][i]
        pp = WWZ_simple_linear[2][i]
        de = saca_punto(ff, pp, 6)
    
        l1.append(de[0])
    
        if len(de) >= 2:
            l2.append(de[1])
        else:
            l2.append(None) 
        
        if len(de) >= 3:
            l3.append(de[2])
        else:
            l3.append(None)  
        
        if len(de) >= 4:
            l4.append(de[3])
        else:
            l4.append(None)  

        if len(de) >= 5:
            l5.append(de[4])
        else:
            l5.append(None)
    return l1,l2,l3,l4,l5


def alaWWZ2(WWZ_simple_linear):
    empty_matrices = [[] for _ in range(len(WWZ_simple_linear[0]))]
    
    for i in range(len(WWZ_simple_linear[0])):
        ff = WWZ_simple_linear[1][i]
        pp = WWZ_simple_linear[2][i]
        dede = saca_punto2(ff, pp)
        
        for row in dede:
            row.append(WWZ_simple_linear[0][i][0])
            
        empty_matrices[i] = dede

    return empty_matrices


def puntos_b(WWZ_simple_linear):
    matrices_list = alaWWZ2(WWZ_simple_linear)
    big_matrix = np.vstack(matrices_list)
    return big_matrix


############################################################################################################################################################################################################################

def wwz_avep_p(WWZ_simple_linear,xx=0):
    df = pd.DataFrame(WWZ_simple_linear[2])
    df2 = pd.DataFrame(WWZ_simple_linear[1])
    df11=df.mean()
    df22=df2.mean()
    freq1 = df22.to_numpy() 
    Averag = df11.to_numpy() 

    plt.figure(figsize=(15,6))
    plt.plot(1/freq1, Averag, 'ro',markersize=3,color='r')
    plt.plot(1/freq1, Averag, '-',markersize=3,alpha=0.4,color='r')
    # line of search
    plt.axvline(x=xx, linewidth=2, color='b')
    plt.xlabel('Period')
    plt.ylabel('Average Power')
    plt.axis([np.min(1/freq1), np.max(1/freq1), 0, np.max(Averag)])
    plt.show()

    return freq1,Averag

def wwz_avep_f(WWZ_simple_linear,xx=0):
    df = pd.DataFrame(WWZ_simple_linear[2])
    df2 = pd.DataFrame(WWZ_simple_linear[1])
    df11=df.mean()
    df22=df2.mean()
    freq1 = df22.to_numpy() 
    Averag = df11.to_numpy() 

    plt.figure(figsize=(15,6))
    plt.plot(freq1, Averag, 'ro',markersize=3,color='r')
    plt.plot(freq1, Averag, '-',markersize=3,alpha=0.4,color='r')
    # line of search
    plt.axvline(x=xx, linewidth=2, color='b')
    plt.xlabel('Freq')
    plt.ylabel('Average Power')
    plt.axis([np.min(freq1), np.max(freq1), 0, np.max(Averag)])
    plt.show()

    return freq1,Averag

def wwz_avep_d_f(WWZ_simple_linear,xx=0):
    df = pd.DataFrame(WWZ_simple_linear[2])
    df2 = pd.DataFrame(WWZ_simple_linear[1])
    df11=df.mean()
    df22=df2.mean()
    dd33=df.std()
    dd44=df.var()
    
    
    freq1 = df22.to_numpy() 
    Averag = df11.to_numpy() 
    desvi = dd44.to_numpy() 
    
    plt.figure(figsize=(15,6))

    # line of search
    #plt.plot(freq1, Averag, 'ro',markersize=3,color='r')
    
    plt.plot(freq1, Averag, '-',markersize=3,alpha=0.4,color='r')
    plt.scatter(freq1, Averag, c=desvi, cmap='coolwarm', s=50, edgecolors=None)
    # line of search
    plt.axvline(x=xx, linewidth=2, color='b')
    plt.xlabel('Freq')
    plt.ylabel('Average Power')
    plt.axis([np.min(freq1), np.max(freq1), 0, np.max(Averag)])
    plt.show()

    return freq1,Averag,desvi
############################################################################################################################################################################################################################

def lomb_sc(time,signal,mini=0.01,maxi=10,step=0.001):

    signal2 =  signal-signal.mean()

    # GLS with frequency 0.01 to 10 and 0.001 step
    rhool = LombScargle(time, signal)
    frequency = np.arange(mini, maxi, step)
    power = rhool.power(frequency)

    # obtain max (peak) frequency
    index_of_peak = np.argmax(power)
    period_of_peak = 1 / frequency[index_of_peak]

    #residual (signal - model of frequency peak)
    residual = signal - LombScargle(time, signal).model(time, frequency[index_of_peak])


    phase = (time / period_of_peak) % 1


    plt.figure(figsize=(10, 8))
    plt.subplot(4, 1, 1)
    plt.scatter(time, signal, label='Original', s=5)
    plt.xlabel('Time')
    plt.ylabel('Signal')
    plt.legend()


    plt.subplot(4, 1, 2)
    plt.plot(frequency, power, label='Lomb Scargle')
    plt.xlim(0,maxi)
    plt.xlabel('Freq')
    plt.ylabel('Power')
    plt.legend()

    plt.subplot(4, 1, 3)
    plt.scatter(time, residual, label='Residual', s=5)
    plt.xlabel('Time')
    plt.ylabel('Residual')
    plt.legend()

    con_phase = pd.concat([phase, phase+1])
    con_signal = pd.concat([signal, signal])
    con_sin = np.concatenate((LombScargle(time, signal).model(time, frequency[index_of_peak]), LombScargle(time, signal).model(time, frequency[index_of_peak])))

    plt.subplot(4, 1, 4)
    plt.scatter(con_phase, con_signal, s=5)
    plt.scatter(con_phase, con_sin, s=5)
    plt.xlabel('Phase')
    plt.ylabel('Signal')


    plt.tight_layout()
    plt.show()

    print(f"PEAK of FREQ {frequency[index_of_peak]}")
    
    return residual,frequency[index_of_peak]

def lomb_sc2(time,signal,mini=0.01,maxi=10,step=0.001):

    signal2 =  signal-signal.mean()
    # GLS with frequency 0.01 to 10 and 0.001 step
    rhool = LombScargle(time, signal)
    frequency = np.arange(mini, maxi, step)
    power = rhool.power(frequency)

    # obtain max (peak) frequency
    index_of_peak = np.argmax(power)
    period_of_peak = 1 / frequency[index_of_peak]

    #residual (signal - model of frequency peak)
    residual = signal - LombScargle(time, signal).model(time, frequency[index_of_peak])


    phase = (time / period_of_peak) % 1


    plt.figure(figsize=(10, 8))
    plt.subplot(4, 1, 1)
    plt.scatter(time, signal, label='Original', s=5)
    plt.xlabel('Time')
    plt.ylabel('Signal')
    plt.legend()


    plt.subplot(4, 1, 2)
    plt.plot(frequency, power, label='Lomb Scargle')
    plt.xlim(0,maxi)
    plt.xlabel('Freq')
    plt.ylabel('Power')
    plt.legend()

    plt.subplot(4, 1, 3)
    plt.scatter(time, residual, label='Residual', s=5)
    plt.xlabel('Time')
    plt.ylabel('Residual')
    plt.legend()

    con_phase = pd.concat([phase, phase+1])
    con_signal = pd.concat([signal, signal])
    con_sin = np.concatenate((LombScargle(time, signal).model(time, frequency[index_of_peak]), LombScargle(time, signal).model(time, frequency[index_of_peak])))

    plt.subplot(4, 1, 4)
    plt.scatter(con_phase, con_signal, s=5)
    plt.scatter(con_phase, con_sin, s=5)
    plt.xlabel('Phase')
    plt.ylabel('Signal')


    plt.tight_layout()
    plt.show()

    print(f"PEAK of FREQ {frequency[index_of_peak]}")
    
    return residual,frequency[index_of_peak],frequency,power

def lomb_sc_short(time,signal,mini=0.01,maxi=10,step=0.001):

    signal2 =  signal-signal.mean()

    # GLS with frequency 0.01 to 10 and 0.001 step
    rhool = LombScargle(time, signal)
    frequency = np.arange(mini, maxi, step)
    power = rhool.power(frequency)

    # obtain max (peak) frequency
    index_of_peak = np.argmax(power)
    period_of_peak = 1 / frequency[index_of_peak]

    #residual (signal - model of frequency peak)
    residual = signal - LombScargle(time, signal).model(time, frequency[index_of_peak])


    phase = (time / period_of_peak) % 1

    ddd = lista_peak_fp(power, frequency)
    
    plt.figure(figsize=(15,6))
    plt.grid()
    plt.plot(frequency, power, label='Lomb Scargle',c='r')
    plt.scatter(ddd['Freq'],ddd['Power'],marker = '*',c='b')
    plt.xlim(mini,maxi)
    plt.xlabel('Freq')
    plt.ylabel('Power')
    plt.legend() 
    
    
    
    return ddd




def lomb_sc_white_list(time,signal,ite=1,mini=0.01,maxi=10,step=0.001):
    residual = signal
    feq = []
    per = []
    powe = []
    
    for i in range(ite):
        
        signal2 =  residual-residual.mean()

        # GLS with frequency 0.01 to 10 and 0.001 step
        rhool = LombScargle(time, residual,normalization='psd')
        frequency = np.arange(mini, maxi, step)
        power = rhool.power(frequency)

        # obtain max (peak) frequency
        index_of_peak = np.argmax(power)
        period_of_peak = 1 / frequency[index_of_peak]

        residual = residual - LombScargle(time, residual,normalization='psd').model(time, frequency[index_of_peak])


        phase = (time / period_of_peak) % 1

        feq.append(frequency[index_of_peak])
        per.append(1/frequency[index_of_peak])
        powe.append(power[index_of_peak])
    
    resp = pd.DataFrame({'Freq': feq, 'Period': per, 'Power': powe})
    
    return resp


def lomb_normal(time,datos,mini=0.001,maxi=10,step=0.0001):
    rhool = LombScargle(time, datos,normalization='psd')
    frequency = np.arange(mini, maxi, step)
    power = rhool.power(frequency)
    
    return frequency,power








############################################################################################################################################################################################################################


def PDM_B(x,y,mini=0.01,maxi=1,step=0.01,tipo='period'):
    
    S = pyPDM.Scanner(minVal=mini, maxVal=maxi, dVal=step, mode=tipo)

    P = pyPDM.PyPDM(x, y)

    f1, t1 = P.pdmEquiBinCover(10, 3, S)

    plt.figure(facecolor='white')
    plt.title("Result of PDM analysis")
    plt.xlabel("Frequency")
    plt.ylabel("Theta")
    plt.plot(f1, t1, 'bp-')
    plt.legend(["pdmEquiBinCover", "pdmEquiBin"])
    plt.show()
    posicion_minimo = np.argmin(t1)
    valor_correspondiente = f1[posicion_minimo]
    
    
    
    return valor_correspondiente

def PDM_G(x,y,mini=0.01,maxi=1,step=0.01,tipo='period'):
    
    S = pyPDM.Scanner(minVal=mini, maxVal=maxi, dVal=step, mode=tipo)

    P = pyPDM.PyPDM(x, y)

    f1, t1 = P.pdmEquiBin(10, S)

    plt.figure(facecolor='white')
    plt.title("Result of PDM analysis")
    plt.xlabel("Frequency")
    plt.ylabel("Theta")
    plt.plot(f1, t1, 'bp-',c='g')
    plt.legend(["pdmEquiBinCover", "pdmEquiBin"])
    plt.show()
    posicion_minimo = np.argmin(t1)
    valor_correspondiente = f1[posicion_minimo]
    
    
    
    return valor_correspondiente





############################################################################################################################################################################################################################
############################################################################################################################################################################################################################

def wwz_beta(tt,dd,N=100,freq_low=0.1,freq_high=2,freq_steps=0.01,b=1):
    
    f = 2
    w = 2 * np.pi * f
    c = 1/(2*w**2)
    freq_lin = [freq_low, freq_high, freq_steps]
    
    wwz_result = wwz.wwt(tt, dd, N, freq_lin, b*c, 'linear')
    
    return wwz_result


def wwz_plot(wwz_result):
    pnt = puntos_b(wwz_result)
    inten = pnt[:, 1]

    fig, ax = plt.subplots(figsize=(15,6))
    heatmap = ax.pcolormesh(wwz_result[0], wwz_result[1], wwz_result[2],cmap='RdYlGn')


    ax.set_xlabel('Time')
    ax.set_ylabel('Freq')

    cbar = plt.colorbar(heatmap)

    plt.scatter(pnt[:, 2],pnt[:, 0],color='black',s=4,alpha=pnt[:, 1]/np.max(inten))

    plt.show()
    
    return

def wwz_beta2(tt,dd,N=100,freq_low=0.1,freq_high=2,freq_steps=0.01,b=1):
    
    f = 2
    w = 2 * np.pi * f
    c = 1/(2*w**2)
    freq_lin = [freq_low, freq_high, freq_steps]
    
    wwz_result = wwz.wwt(tt, dd, N, freq_lin, b*c, 'periodo')
    
    return wwz_result


def wwz_plot2(wwz_result):

    fig, ax = plt.subplots(figsize=(15,6))
    heatmap = ax.pcolormesh(wwz_result[0], wwz_result[1], wwz_result[2],cmap='magma')


    ax.set_xlabel('Time')
    ax.set_ylabel('Freq')

    cbar = plt.colorbar(heatmap)


    plt.show()
    
    return

def wwz_fast(tt,dd,N=100,freq_low=0.1,freq_high=2,freq_steps=0.01,b=1):
    
    wwz = WWZ(tt, dd, 0)

    freq = wwz.get_freq(p_min=1/freq_high, p_max=1/freq_low,diff=freq_steps)
    tau = wwz.get_tau(n_bins=N)
    wwz.set_freq(freq,verbose=False)
    wwz.set_tau(tau,verbose=False)    
    wwz.transform(c=b*0.012665147955292222,verbose=1)
    pass
    valores = np.arange(freq_low+freq_steps, freq_high+freq_steps, freq_steps)
    numero_de_copias = N
    freq_apilado = np.tile(valores,  (N, 1))
    la = np.shape(freq_apilado)
    timess = np.linspace(tt[0], tt[-1], N)
    vectores = []
    
    for t in timess:
        vector = np.full(la[1], t)
        vectores.append(vector)

    time_apilado = np.vstack(vectores)
    
    
    result_wwz = wwz.wwz
    result_wwa = wwz.wwa
    
    return [time_apilado,freq_apilado,result_wwz,result_wwa]

def wwz_fast2(tt,dd,N=100,freq_low=0.1,freq_high=2,freq_steps=0.01,b=0.25):
    
    wwz = WWZ(tt, dd, 0)

    freq = wwz.get_freq(p_min=1/freq_high, p_max=1/freq_low,diff=freq_steps)
    tau = wwz.get_tau(n_bins=N)
    wwz.set_freq(freq,verbose=False)
    wwz.set_tau(tau,verbose=False)    
    wwz.transform(c=b*0.012665147955292222,verbose=1)
    pass
    valores = np.arange(freq_low, freq_high, freq_steps)
    numero_de_copias = N
    freq_apilado = np.tile(valores,  (N, 1))
    la = np.shape(freq_apilado)
    timess = np.linspace(tt[0], tt[-1], N)
    vectores = []
    
    for t in timess:
        vector = np.full(la[1], t)
        vectores.append(vector)

    time_apilado = np.vstack(vectores)
    
    
    result_wwz = wwz.wwz
    result_wwa = wwz.wwa
    
    return [time_apilado,freq_apilado,result_wwz,result_wwa]
############################################################################################################################################################################################################################
############################################################################################################################################################################################################################

def plot_phase(time,signal,stepp=0.001,val=9,ss=10):
    
    def phaseee(angle_degrees):
        phase = (time / (angle_degrees)) % 1

        con_phase = pd.concat([phase, phase+1])
        con_signal = pd.concat([signal, signal])

        plt.figure(figsize=(12, 6)) 
        plt.scatter(con_phase, con_signal, s=ss,c='m')
        plt.xlabel('Phase')
        plt.ylabel('Signal')
        plt.show()


    interact(phaseee, angle_degrees=BoundedFloatText(min=0, max=1000,step=stepp, value=val));
    


def plot_phase2(time,signal,stepp=0.001,val=9,ss=10):
    
    aaa = signal.max()
    bbb = signal.min()
    
    ccc = 0.5 * (aaa - bbb)
    
    aa = aaa + ccc
    bb = bbb - ccc
    
    def phaseee(angle_degrees):
        phase = (time / (angle_degrees)) % 1

        con_phase = pd.concat([phase, phase+1])
        con_signal = pd.concat([signal, signal])

        plt.figure(figsize=(12, 6)) 
        plt.scatter(con_phase, con_signal, s=ss,c='m')
        plt.ylim(aa,bb)
        plt.xlabel('Phase')
        plt.ylabel('Signal')
        plt.show()


    interact(phaseee, angle_degrees=BoundedFloatText(min=0, max=1000,step=stepp, value=val));


    
    
################################################################################################################################################################################################################################################



def ruido_rojo(time,lt,mini=0.0001,maxi=100,step=0.0001,methodss='lm'):
    
    def funci_log(f, a, b,c,d):
        return (d/(1+ (2*np.pi*c*f)**a)) + b 

    lombs = lomb_normal(time,lt,mini=0.0001,maxi=100,step=0.0001)
    
    parametros, covariance35_39 = curve_fit(funci_log, lombs[0], np.sqrt(np.abs(lombs[1]))* np.sqrt(4.0 / len(time)),method=methodss)
    
    plt.figure(figsize=(24, 8)) 
    plt.scatter(lombs[0], np.sqrt(np.abs(lombs[1]))* np.sqrt(4.0 / len(time)), label='Lomb Scargle', color='r',s=3)
    plt.plot(lombs[0], funci_log(lombs[0], *parametros), 'g-', label='curve')
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel('Log freq')
    plt.ylabel('Log AMP')
    plt.legend()
    plt.show()
    argu = np.array(['Gamma', 'Alpha 0', 'Tau', 'Alpha w'])
    
    df = pd.DataFrame({'L': argu,'Param': parametros})
    df
    
    return df,lombs,[lombs[0], funci_log(lombs[0], *parametros)]




################################################################################################################################################################################################################################################



def ruido_rojo2(time,lt,mini=0.0001,maxi=100,step=0.0001,methodss='lm'):
    
    def funci_log(f, a, b,c,d):
        return (d/(1+ (2*np.pi*c*f)**a)) + b 

    lombs = lomb_normal(time,lt,mini,maxi,step)
    
    parametros, covariance35_39 = curve_fit(funci_log, lombs[0], np.sqrt(np.abs(lombs[1]))* np.sqrt(4.0 / len(time)),method=methodss)
    
    plt.figure(figsize=(24, 8)) 
    plt.scatter(lombs[0], np.sqrt(np.abs(lombs[1]))* np.sqrt(4.0 / len(time)), label='Lomb Scargle', color='r',s=3)
    plt.plot(lombs[0], funci_log(lombs[0], *parametros), 'g-', label='curve')
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel('Log freq')
    plt.ylabel('Log AMP')
    plt.legend()
    plt.show()
    argu = np.array(['Gamma', 'Alpha 0', 'Tau', 'Alpha w'])
    
    df = pd.DataFrame({'L': argu,'Param': parametros})
    df
    
    return df,lombs,[lombs[0], funci_log(lombs[0], *parametros)]



################################################################################################################################################################################################################################################



def analizar_frecuencias_ald(frequencies, max_error,output_file="frecuencia.dat"):
    related_frequencies = set()
    relations_dict = {}

    def verificar_relaciones(operation):
        nonlocal relations_dict
        if operation == 'multiplo':
            print("\nRelaciones de múltiplo:")
            for i, freq1 in enumerate(frequencies):
                if i in related_frequencies:
                    continue
                for j in range(i):
                    freq2 = frequencies[j]
                    multiple = freq1 / freq2
                    error = abs(freq1 - round(multiple) * freq2)
                    if error <= max_error:
                        print(f"F{i} ~ {round(multiple)} F{j}")
                        related_frequencies.add(i)
                        relation_key = f"F{i}"
                        relation_value = f"{round(multiple)} F{j}"
                        if relation_key in relations_dict:
                            relations_dict[relation_key].append(relation_value)
                        else:
                            relations_dict[relation_key] = [relation_value]
                        break

        elif operation == 'suma' or operation == 'resta':
            print(f"\nRelaciones de {operation}:")
            for k in range(2, len(frequencies)):  # bucle de k iniciando en k=2
                if k in related_frequencies:
                    continue
                freq3 = frequencies[k]
                best_relation = None
                best_error = float('inf')
                for i in range(k):
                    if i in related_frequencies:
                        continue
                    freq1 = frequencies[i]
                    for j in range(k):
                        if j in related_frequencies or i == j:  # Ignorar cuando i = j
                            continue
                        freq2 = frequencies[j]
                        for a in range(1, 10):
                            for b in range(1, 10):
                                if operation == 'suma':
                                    error = abs(freq3 - (freq1 * a + freq2 * b))
                                else:  # operation == 'resta'
                                    result = freq1 * a - freq2 * b
                                    if result >= 0:
                                        error = abs(freq3 - result)
                                    else:
                                        # error = abs(freq3 + result) 
                                        error = float('inf') 
                                if error <= max_error and error < best_error:
                                    best_error = error
                                    best_relation = (a, b, i, j)  # Almacenamos también los índices de las frecuencias i y j
                if best_relation:
                    a, b, i, j = best_relation  # Desempaquetamos los valores de la mejor relación
                    if (operation == 'suma' and a > 0 and b > 0) or \
                            (operation == 'resta' and b > 0 and a > 0):
                        print(f"F{k} ~  {a}F{i} {'+' if operation == 'suma' else '-'} {b}F{j}")
                        related_frequencies.add(k)
                        relation_key = f"F{k}"
                        relation_value = f"{a}F{i} {'+' if operation == 'suma' else '-'} {b}F{j}"
                        if relation_key in relations_dict:
                            relations_dict[relation_key].append(relation_value)
                        else:
                            relations_dict[relation_key] = [relation_value]

    verificar_relaciones('multiplo')
    verificar_relaciones('suma')
    verificar_relaciones('resta')

    print("\nFrecuencias independientes:")
    for idx, freq in enumerate(frequencies):
        if idx not in related_frequencies:
            print(f"F{idx} es independiente.")

    # Escribir relaciones en un archivo
    with open(output_file, "w") as file:
        file.write("Frecuencia\tRelación\n")
        for key, value in sorted(relations_dict.items(), key=lambda item: int(item[0][1:])):  # Ordenar por nombre de F
             for relation in value:
                file.write(f"{key}\t{relation}\n")

        # Escribir frecuencias independientes
        for idx, freq in enumerate(frequencies):
            if idx not in related_frequencies:
                file.write(f"F{idx}\tIndep\n")




































