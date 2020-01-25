import numpy as np
import matplotlib.pyplot as plt



def plotPSD(xs=[],ys=[]):
    dts = []
    for i in range(1, len(xs), 1):
        dts.append(xs[i] - xs[i - 1])
    dtAvg = np.average(dts)
    dtStdDev = np.std(dts)
    fSamp = 1.0/dtAvg
    sig = ys
    s = sig



    fft = np.fft.fft(s)
    T = dtAvg  # sampling interval
    N = len(s)

    # 1/T = frequency
    f = np.linspace(0, 1 / T, N)

    t1 = f[:N // 2]
    t2 = np.abs(fft)[:N // 2] * 1 / N

    for i in range(0,len(t1),1):
        print('%s        %s'%(t1[i],t2[i]))

    plt.clf()
    plt.figure(1)
    plt.plot(t1,t2)
    plt.show()
