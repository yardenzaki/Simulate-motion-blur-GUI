##########################################################################
# program: psd_syn.py
# author: Tom Irvine
# version: 1.5
# date: September 13, 2013
# description:  this script synthesis a time history to satisfy a
#               power spectral density
#
##########################################################################

from __future__ import print_function

from toms_py import signal_stats, sample_rate_check, small
from toms_py import enter_float
from toms_py import GetInteger2
from toms_py import time_history_plot
from toms_py import WriteData2

from signal_utilities import EnterPSD

from velox_correction import velox_correction

from math import pi, ceil

import numpy
from numpy import zeros, random, linspace, std
from numpy import sqrt, log, cos, sin, log10
from numpy import argmax, interp

from scipy.fftpack import fft
from scipy.fftpack import ifft
from scipy.stats import norm

import matplotlib.pyplot as plt

from sys import stdin


########################################################################

def magnitude_resolve(mmm, mH, Y):
    #
    mHm1 = mH - 1
    z = zeros(mH, 'f')
    mag_seg = zeros(mH, 'f')
    #
    #     for i in range (0,mH):
    #       z[i]=sqrt(Y.real[i]**2+Y.imag[i]**2)
    #
    z = abs(Y) / float(mmm)
    #
    mag_seg[0] = z[0] ** 2
    #
    mag_seg[1:mHm1] = ((2 * z[1:mHm1]) ** 2) / 2
    #
    return mag_seg


########################################################################

def Hanning_initial(mmm):
    H = zeros(mmm, 'f')
    tpi = 2 * pi
    alpha = linspace(0, tpi, mmm)
    ae = sqrt(8. / 3.)
    H = ae * 0.5 * (1. - cos(alpha))
    return H


########################################################################

def psd_core(mmm, mH, maxf, NW, b,df,H):
    full = zeros(mH, 'f')

    mag_seg = zeros(mH, 'f')
    amp_seg = zeros(mmm, 'f')
    den = df * (2 * NW - 1)

    nov = 0
    for ijk in range(1, 2 * NW):
        amp_seg[0:mmm] = b[(0 + nov):(mmm + nov)]

        nov = nov + int(mmm / 2)

        mean = sum(amp_seg) / float(mmm)
        amp_seg -= mean

        amp_seg *= H

        Y = fft(amp_seg)

        mag_seg = magnitude_resolve(mmm, mH, Y)

        full += mag_seg

    full /= den

    ms = sum(full)

    return full, ms


########################################################################

def psd_syn_main():
    tpi = 2 * pi

    freq_spec, amp_spec, rms, num, slope, FileName = EnterPSD()

    nm1 = num - 1
    LS = nm1

    three_rms = 3 * rms

    print(" ")
    print(" Enter duration(sec)")
    tmax = enter_float()

    fmax = max(freq_spec)

    #Control the sample rate here
    sample_rate_factor=10 #so that the sample rate will be 10 times the maximal freq in the PSD
    #sr = fmax * 20. #OLD
    sr = fmax * sample_rate_factor
    dt = 1 / sr

    spec_grms = rms

    np = ceil(tmax / dt)
    np3 = 3 * np
    print("np3", np3)
    mu = 0
    sigma = 1

    print(" Generate white noise")

    white_noise = random.normal(mu, sigma, np3)

    print(" End white noise")

    num_fft = 2

    while (num_fft < np):
        num_fft *= 2

    N = num_fft
    df = 1. / (N * dt)

    m2 = int(num_fft / 2)

    fft_freq = linspace(0, (m2 - 1) * df, m2)
    fft_freq2 = linspace(0, (num_fft - 1) * df, num_fft)

    spec = zeros(m2, float)
    sq_spec = zeros(m2, float)

    js = 0

    print(" Interpolate specification")

    if (fft_freq[0] <= 0):
        fft_freq[0] = 0.5 * fft_freq[1]

    x = log10(fft_freq)
    xp = log10(freq_spec)
    yp = log10(amp_spec)

    y = interp(x, xp, yp, left=-10, right=-10)

    sq_spec = sqrt(10 ** y)

    # add option for sine tones later

    print(" Calculating FFT ")

    white_noise_trunc = white_noise[0:num_fft]

    Y = zeros(num_fft, complex)

    YF = fft(white_noise_trunc)

    print(" Apply spec ")

    # for j in range(1, m2):
    #    Y[j] = sq_spec[j] * YF[j]

    YFn = YF[0:m2]

    Y[0:m2] = sq_spec * YFn

    Y[0] = 0.

    print(" Make symmetric")

    for j in range(1, m2):
        Y[num_fft - j] = complex(Y[j].real, -Y[j].imag)

    print(" Calculating inverse FFT ")

    YI = ifft(Y)

    print(" YIR")

    YIR = YI.real

    print(" psd_th_m2")

    psd_th_m2 = YIR

    nL = m2

    print(" psd_th")

    psd_th = YIR[0:np]

    np = len(psd_th)

    TT = linspace(0, (np - 1) * dt, np)

    print(" ")
    print("num_fft=%d " % num_fft)
    print("np=%d      " % np)

    stddev = std(psd_th)

    psd_th *= (spec_grms / stddev)

    #### check psd ############################################################

    a = TT
    b = psd_th

    num = len(a)

    sr, dt, mean, sd, rms, skew, kurtosis, dur = signal_stats(a, b, num)

    sr, dt = sample_rate_check(a, b, num, sr, dt)

    ########################################################################

    # print " "
    # print " Remove mean:  1=yes  2=no "

    mr_choice = 1

    # print " "
    # print " Select Window: 1=Rectangular 2=Hanning "

    # h_choice = GetInteger2()

    h_choice = 2

    ########################################################################

    n = num

    ss = zeros(n)
    seg = zeros(n, 'f')
    i_seg = zeros(n)

    NC = 0
    for i in range(0, 1000):
        nmp = 2 ** (i - 1)
        if (nmp <= n):
            ss[i] = 2 ** (i - 1)
            seg[i] = float(n) / float(ss[i])
            i_seg[i] = int(seg[i])
            NC = NC + 1
        else:
            break

    print(' ')
    print(' Number of   Samples per   Time per        df    ')
    print(' Segments     Segment      Segment(sec)   (Hz)   dof')

    for i in range(1, NC + 1):
        j = NC + 1 - i
        if j > 0:
            if (i_seg[j] > 0):
                tseg = dt * ss[j]
                ddf = 1. / tseg
                print('%8d \t %8d \t %10.3g  %10.3g    %d' \
                      % (i_seg[j], ss[j], tseg, ddf, 2 * i_seg[j]))
        if (i == 12):
            break

    ijk = 0
    while ijk == 0:
        print(' ')
        print(' Choose the Number of Segments:  ')
        s = stdin.readline()
        NW = int(s)
        for j in range(0, len(i_seg)):
            if NW == i_seg[j]:
                ijk = 1
                break

    # check

    mmm = 2 ** int(log(float(n) / float(NW)) / log(2))

    df = 1. / (mmm * dt)

    # begin overlap

    mH = int(((mmm / 2) - 1))

    print(" ")
    print("     number of segments   NW= %d " % NW)
    print("       samples/segments  mmm= %d " % mmm)
    print(" half samples/segment-1   mH=%d  " % mH)
    print(" ")
    print("        df=%6.3f Hz" % df)

    maxf = (mH - 1) * df

    H = Hanning_initial(mmm)

    freq = zeros(mH, 'f')
    freq = linspace(0, maxf, mH)

    delta = freq_spec[0] / 30

    ########################################################################

    nnt = 3

    print(" Velocity correction")

    for kvn in range(0, nnt):

        acc, velox, dispx = velox_correction(psd_th, dt, freq_spec[0])

        ratio = spec_grms / std(psd_th)

        psd_th *= ratio
        velox *= ratio
        dispx *= ratio

        full, ms = psd_core(mmm, mH, maxf, NW, psd_th,df,H)

        MK = len(psd_th)
        tim = linspace(0, MK * dt, MK)

        # low frequency correction

        psd1 = 0

        for i in range(0, len(freq)):
            #       print " %8.4g  %8.4g " %(freq[i],freq_spec[0])
            if (freq[i] <= freq_spec[0] and freq_spec[0] <= freq[i + 1]):
                x = freq_spec[0] - freq[i]
                c2 = x / (freq[i + 1] - freq[i])
                c1 = 1 - c2
                psd1 = c1 * full[i] + c2 * full[i + 1]
                break

        # print ("\n @@ kvn=%d  psd1=%8.4g  amp_spec=%8.4g  " %(kvn,psd1,amp_spec[0]))

        if (psd1 < amp_spec[0]):
            ca = sqrt(2) * sqrt(amp_spec[0] * df - psd1 * df)

            print("kvn=%d    ca=%9.5g " % (kvn, ca))

            pha = tpi * random.random()

            fff = freq_spec[0] + (-0.5 + random.random()) * delta

            if (kvn == 0):
                fff = freq_spec[0]

            if (kvn == 1):
                fff = freq_spec[0] + delta / 2

            if (kvn == 2):
                fff = freq_spec[0] - delta / 2

            psd_th += ca * sin(tpi * fff * tim + pha)
        ##factoring Time History values to satisfy the PSD factoring
        #We need to factor because:
        #1. we factored LOW PSD to HIGH PSD by multiplying by 1e12
        #2. we synthesized TH based on the HIGH PSD so we have HIGH TH
        #3. Now we factor the HIGH TH to LOW TH by dividing by 1e6 (sqrt of PSD factor)
        psd_th=psd_th*1e-6 #now we are in the original units which are rad
        # print("psd_th",psd_th)
        # print("type", type(psd_th))
    ################################################################################

    tempf = freq[0:mH - 1]
    tempa = full[0:mH - 1]
    freq = tempf
    full = tempa

    rms = sqrt(ms * df)
    three_rms = 3 * rms

    psd_rms = rms

    print(" ")
    print(" Overall RMS = %10.3g " % rms)
    print(" Three Sigma = %10.3g " % three_rms)

    idx = argmax(full)

    print(" ")
    print(" Maximum:  Freq=%8.4g Hz   Amp=%8.4g " % (freq[idx], full[idx]))

    ########################################################################

    print(" ")
    print(" ** Output filenames ** ")

    print(" ")
    print("Writing time history .... ")
    #a_output_file_path = stdin.readline()
    a_output_file_path="Out Time History "+ FileName +" Rad"+".txt"

    # print(" ")
    # print("Enter the velocity time history filename ")
    # v_output_file_path = stdin.readline()
    #
    # print(" ")
    # print("Enter the displacement time history filename ")
    # d_output_file_path = stdin.readline()

    file_path = a_output_file_path.rstrip('\n')
    WriteData2(np, TT, psd_th, file_path)

    # file_path = v_output_file_path.rstrip('\n')
    # WriteData2(np, TT, velox, file_path)
    #
    # file_path = d_output_file_path.rstrip('\n')
    # WriteData2(np, TT, dispx, file_path)

    print(" ")
    print(" Writing PSD data to file")
    #iacc = GetInteger2()
    iacc=1

    if (iacc == 1):
        print(" ")
        print("Writing output PSD ... ")
        output_file_path = "Out PSD " + FileName +".txt"
        output_file = output_file_path.rstrip('\n')
        mH = len(freq)
        WriteData2(mH, freq, full, output_file)

    ###########################################################################

    n = len(TT)

    max_lines = 400000

    psd_th_store = psd_th

    if (n > max_lines):
        print(' ')
        K = int(ceil(n / max_lines))
        TTo = TT
        print('processing acceleration plot data')
        TT, psd_th = small(TTo, psd_th, K)
        # print('processing velocity plot data')
        # TT, velox = small(TTo, velox, K)
        # print('processing displacement plot data')
        # TT, dispx = small(TTo, dispx, K)

    string_value1 = "  %6.3g" % sd

    title_string = FileName+' Time History based on PSD' + string_value1 + 'X1e-6 RMS'
    time_history_plot(TT, psd_th, 1, 'Time(sec)', 'Angle(rad)', title_string, FileName+' Time History')

    # sd = std(velox)
    # string_value1 = "  %6.3g" % sd
    #
    # title_string = 'Velocity' + string_value1 + ' in/sec RMS'
    # time_history_plot(TT, velox, 2, 'Time(sec)', 'Vel(in/sec)', title_string, 'velox')
    #
    # sd = std(dispx)
    # string_value1 = "  %6.3g" % sd
    #
    # title_string = 'Displacement' + string_value1 + ' inch RMS'
    # time_history_plot(TT, dispx, 3, 'Time(sec)', 'Disp(inch)', title_string, 'disp')

    ############################################################################

    plt.figure(4)
    psd_th_urad = psd_th * 1e6 #switch rad to urad
    #print("psd_th_urad", psd_th_urad, type(psd_th_urad))
    maxVal=psd_th_urad.max()
    minVal=psd_th_urad.min()
    mu0, std0 = norm.fit(psd_th_urad)
    #print("maxVal", maxVal, minVal)
    plt.hist(psd_th_urad, bins=51,density=True)
    xmin, xmax = plt.xlim()
    xx = numpy.linspace(xmin, xmax, 100)
    pp = norm.pdf(xx, mu0, std0)
    plt.plot(xx, pp, 'k', linewidth=2)


    plt.xlim(xmin=-maxVal, xmax=maxVal)
    plt.ylim(ymin=0, ymax=0.1)
    #plt.yscale('log')
    title0= FileName+" Tilts Histogram; Fit Values: {:.2f} and {:.2f}".format(mu0, std0)
    plt.title(title0)
    plt.xlabel(' Angel (urad)')
    plt.ylabel(' Counts ')
    plt.savefig(FileName + ' Hist')
    ############################################################################

    pmin = 10 ** 40
    pmax = 10 ** -40

    fmin = 10 ** 40
    fmax = 10 ** -40

    fmin = min(freq_spec)
    fmax = max(freq_spec)

    for i in range(0, len(freq)):
        if freq[i] > fmax:
            break
        if full[i] > 0 and freq[i] >= fmin and full[i] > pmax:
            pmax = full[i]
        if full[i] > 0 and freq[i] >= fmin and full[i] < pmin:
            pmin = full[i]

    xmax = 10 ** -30
    xmin = xmax

    for i in range(-30, 30):
        if (fmax < 10 ** i):
            xmax = 10 ** i
            break

    for i in range(30, -30, -1):
        if (fmin > 10 ** i):
            xmin = 10 ** i
            break

    ymax = 10 ** -30
    ymin = ymax

    for i in range(-30, 30):
        if (pmax < 10 ** i):
            ymax = 10 ** i
            break

    for i in range(30, -30, -1):
        if (pmin > 10 ** i):
            ymin = 10 ** i
            break

    ###############################################################################

    psd_label = ' (rad^2/Hz)'

    plt.gca().set_autoscale_on(False)

    plt.figure(5)
    plt.plot(freq, full, label="synthesis (magnified by 1*E12)")
    plt.plot(freq_spec, amp_spec, label="specification (magnified by 1*E12)")
    title_string = FileName + ' Power Spectral Density ' + str("%6.3g" % psd_rms) + '*1E-6 RMS '
    plt.title(title_string)
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.ylabel(psd_label)
    plt.xlabel(' Frequency (Hz) ')
    plt.grid(True)
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc="upper left")
    plt.savefig(FileName + ' power_spectral_density')

    ############################################################################

    print(" ")
    print(" view plots ")

    plt.draw()
    plt.show()

    #return TT,psd_th_urad #Time in sec and TH in urad
    return TT, psd_th  # Time in sec and TH in rad