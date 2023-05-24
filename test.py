from obspy.core import read
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import mlab
from matplotlib.colors import Normalize

def _nearest_pow_2(x):
    import math
    a = math.pow(2, math.ceil(np.log2(x)))
    b = math.pow(2, math.floor(np.log2(x)))
    if abs(a - x) < abs(b - x):
        return a
    else:
        return b


if __name__ == "__main__":
    st = read("NE68.mseed")
    st1 = read("NE68_normal.mseed")
    print(len(st))
    ind = 25
    tr = st[ind]
    tr1 = st1[ind]
    print(len(tr))
    print(tr.stats)
    print("-----")
    print(tr.stats.sampling_rate)

    #st[0].spectrogram(log=True)
    #st1[0].spectrogram(log=True)

    y = np.fft.fft(tr.data)

    print()
    res = np.fft.ifft(np.log(np.abs(y)))

    ifs = []
    for i in range(len(res)):
        ifs.append(res[i] * res[i]/len(res))
    # plt.specgram(y, Fs=100)
    # plt.show()

    data = tr.data
    samp_rate = 100.0
    samp_rate = float(samp_rate)
    wlen = 128 / samp_rate
    npts = len(data)
    nfft = int(_nearest_pow_2(wlen * samp_rate))
    if npts < nfft:
        msg = (f'Input signal too short ({npts} samples, window length '
               f'{wlen} seconds, nfft {nfft} samples, sampling rate '
               f'{samp_rate} Hz)')
        raise ValueError(msg)
    mult = int(_nearest_pow_2(8.0))
    mult = mult * nfft
    nlap = int(nfft * float(0.9))
    data = data - data.mean()
    end = npts / samp_rate
    specgram, freq, time = mlab.specgram(data, Fs=samp_rate, NFFT=nfft,
                                         pad_to=mult, noverlap=nlap)
    #specgram, freq, time = mlab.specgram(tr.data, Fs=100)
    specgram = np.sqrt(specgram[1:, :])
    freq = freq[1:]
    _range = float(specgram.max() - specgram.min())
    vmin, vmax = [0.0,1.0]
    vmin = specgram.min() + vmin * _range
    vmax = specgram.min() + vmax * _range
    norm = Normalize(vmin, vmax, clip=True)
    halfbin_time = (time[1] - time[0]) / 2.0
    halfbin_freq = (freq[1] - freq[0]) / 2.0
    log = True
    if log:

        freq = np.concatenate((freq, [freq[-1] + 2 * halfbin_freq]))
        time = np.concatenate((time, [time[-1] + 2 * halfbin_time]))
        # center bin
        time -= halfbin_time
        freq -= halfbin_freq
        #plt.specgram(specgram, Fs=100, mode='psd')
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.pcolormesh(time, freq, specgram,norm=norm)
        ax.set_yscale('log')
        plt.show()

    # fig = plt.figure()
    # ax = fig.add_subplot(2, 1, 1)
    # st[0].spectrogram(log=True, axes=ax)
    # #ax.plot(tr.times("matplotlib"), tr.data, "b-")
    
    # ax1 = fig.add_subplot(2, 1, 2)
    # st1[0].spectrogram(log=True, axes=ax1)
    # #ax1.plot(tr1.times("matplotlib"), tr1.data, "b-")
    # #fig.autofmt_xdate()
    # plt.show()