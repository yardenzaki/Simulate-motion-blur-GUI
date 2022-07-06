"""
Creating a simulation video from frame and LoS error analysis results.
The desired ROI will vibrate according to the LoS error and statistical analysis
The analysis neglect any high - order image artifacts (defocus,astigmatism,coma etc..) and simulates only rigid body tilts ( first order Zernike polynoms)
How does it work?
1. Getting tilts LoS PSD Data [rad^2/Hz vs. Hz].
2. Converting to displacements LoS PSD with w.r.t the detector plane by multiplying with the (focal len)^2
3. Converting it to time history equivalent tilts [rad vs. time].
4. Render the corresponding image motion based on the results.

** V01 - 06/07/22 What's new?
    1. Algorithm improvement: Sampling first DX , DY in "exposure time" based on CDF matching and not in a pure random - leads to non-zero correlation between DX , DY
    2. Post Processing: Making animation of PDF , CDF plots of DX,DY for the 50 first frames.

@Author: Yarden Zaki
@Date: 07/01/2022
@Version: 1.0
@Credits: Tom Irvine for the signal processing codes - vibrationdata.com
@Links: https://github.com/yardenzaki
@License: MIT

"""

import cv2
import os
from matplotlib.pylab import *
import psd_syn
from tompy import enter_float, signal_stats, WriteData2

# determine what OpenCV version we are using
try:
    import cv2.cv as cv

    USE_CV2 = True
except ImportError:
    # OpenCV 3.x does not have cv2.cv submodule
    USE_CV2 = False

import sys
import numpy as np
import tkinter as tk
import tkinter.messagebox as tkMessageBox
from scipy.stats import norm
import scipy.stats as sps
from statsmodels.distributions.empirical_distribution import ECDF
import imageio.v2 as imageio


def quit(root):
    root.destroy()


def read_data():
    if not durr.get():  # do something
        tkMessageBox.showinfo("Warning", "Enter duration", parent=button_read)
        return

    self.tmax = float(self.durr.get())

    """
    f = frequency column
    a = PSD column
    num = number of coordinates
    slope = slope between coordinate pairs    
    """

    print(" ")
    print(" The input file must have two columns: freq(Hz) & accel(G^2/Hz)")

    f, a, num = read_two_columns_from_dialog('Select Input File', self.master)

    print("\n samples = %d " % num)

    f = array(f)
    a = array(a)

    self.maxf = max(f)

    nm1 = num - 1

    slope = zeros(nm1, 'f')

    ra = 0

    for i in range(0, int(nm1)):
        #
        s = log(a[i + 1] / a[i]) / log(f[i + 1] / f[i])

        slope[i] = s
        #
        if s < -1.0001 or s > -0.9999:
            ra += (a[i + 1] * f[i + 1] - a[i] * f[i]) / (s + 1.)
        else:
            ra += a[i] * f[i] * log(f[i + 1] / f[i])

    omega = 2 * pi * a

    av = zeros(num, 'f')
    ad = zeros(num, 'f')

    for i in range(0, int(num)):
        av[i] = a[i] / omega[i] ** 2

    rv = 0

    for i in range(0, int(nm1)):
        #
        s = log(av[i + 1] / av[i]) / log(f[i + 1] / f[i])
        #
        if s < -1.0001 or s > -0.9999:
            rv += (av[i + 1] * f[i + 1] - av[i] * f[i]) / (s + 1.)
        else:
            rv += av[i] * f[i] * log(f[i + 1] / f[i])

    for i in range(0, int(num)):
        ad[i] = av[i] / omega[i] ** 2

    rd = 0

    for i in range(0, int(nm1)):
        #
        s = log(ad[i + 1] / ad[i]) / log(f[i + 1] / f[i])
        #
        if s < -1.0001 or s > -0.9999:
            rd += (ad[i + 1] * f[i + 1] - ad[i] * f[i]) / (s + 1.)
        else:
            rd += ad[i] * f[i] * log(f[i + 1] / f[i])

    rms = sqrt(ra)
    three_rms = 3 * rms

    print(" ")
    print(" *** Input PSD *** ")
    print(" ")

    print(" Acceleration ")
    print("   Overall = %10.3g GRMS" % rms)
    print("           = %10.3g 3-sigma" % three_rms)

    self.grms_in = rms
    self.f = f
    self.a = a

    self.rms = rms
    self.freq_spec = f
    self.amp_spec = a
    self.num = num
    self.slope = slope

    self.button_calculate.config(state='normal')

    self.spec_grms = rms

    plt.ion()
    plt.clf()
    plt.figure(self.fig_num)
    self.fig_num += 1
    plt.plot(f, a)
    title_string = 'Power Spectral Density   ' + str("%6.3g" % rms) + ' GRMS Overall '
    plt.title(title_string)
    plt.ylabel(' Accel (G^2/Hz)')
    plt.xlabel(' Frequency (Hz) ')
    plt.grid(which='both')
    plt.savefig('power_spectral_density')
    plt.xscale('log')
    plt.yscale('log')
    # plt.show()
    plt.close("all")

    self.advise(self)
    self.white_noise(self)
    self.calculate_fft(self)
    self.apply_spec(self)
    self.calculate_invfft(self)

    self.button_calculate.config(state='normal')


def Call_load_psd():
    plt.close("all")
    win = tk.Toplevel()
    Load_PSD_Results(win)


def Load_PSD_Results(root):
    global durr

    master = root  # store the parent
    top = tk.Frame(root)  # frame for all class widgets
    top.pack(side='top')  # pack frame in parent's window

    w, h = master.winfo_screenwidth(), master.winfo_screenheight()
    w = int(2. * (w * 0.30))
    h = int(2. * (h * 0.30))
    master.geometry("%dx%d+0+0" % (w, h))

    master.title("PSD to Time History ver 1.0  by Yarden Zaki")
    crow = 0
    hwtext1 = tk.Label(top, text='Synthesis of time history LoS deviation from rotations PSD')
    hwtext1.grid(row=crow, column=0, columnspan=3, padx=8, pady=7, sticky=tk.W)

    crow += 1

    hwtext1 = tk.Label(top, text='The input file must have two columns: Freq [Hz] & LoS [Rad^2/Hz]')
    hwtext1.grid(row=crow, column=0, columnspan=3, padx=8, pady=5, sticky=tk.W)
    crow = crow + 1

    hwtext4 = tk.Label(top, text='Select Output Units')
    hwtext4.grid(row=crow, column=0, columnspan=1, pady=6, sticky=tk.S)

    hwtext_fn = tk.Label(top, text='Enter Duration (sec)')
    hwtext_fn.grid(row=crow, column=1, columnspan=1, padx=14, pady=10, sticky=tk.S)

    crow = crow + 1

    Lb1 = tk.Listbox(top, height=2, exportselection=0)
    Lb1.insert(1, "G, in/sec, in")
    Lb1.insert(2, "G, cm/sec, mm")
    Lb1.grid(row=crow, column=0, pady=2, sticky=tk.N)
    Lb1.select_set(0)

    durr = tk.StringVar()
    durr.set('')
    dur_entry = tk.Entry(top, width=12, textvariable=durr)
    dur_entry.grid(row=crow, column=1, padx=14, pady=1, sticky=tk.N)

    button_read = tk.Button(top, text="Read Input File", command=quit)
    button_read.config(height=3, width=15)
    button_read.grid(row=crow, column=2, columnspan=1, padx=0, pady=2, sticky=tk.N)

    crow = crow + 1

    hwtextadv = tk.Label(top, text='Select Analysis Option')
    hwtextadv.grid(row=crow, column=0, pady=10)

    crow = crow + 1

    myframe = tk.Frame(top)
    myframe.grid(row=crow, column=0, padx=3)
    scrollbar = tk.Scrollbar(myframe)
    scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
    Lba = tk.Listbox(myframe, width=35, yscrollcommand=scrollbar.set)
    Lba.pack()
    scrollbar.config(command=Lba.yview)

    crow = crow + 1

    button_calculate = tk.Button(top, text="Perform Analysis", command=quit)
    button_calculate.config(height=2, width=15, state='disabled')
    button_calculate.grid(row=crow, column=0, columnspan=1, padx=10, pady=10)

    root = master

    button_quit = tk.Button(top, text="Quit", command=lambda root=root: quit(root))
    button_quit.config(height=2, width=15)
    button_quit.grid(row=crow, column=1, padx=6, pady=10)

    button_ex = tk.Button(top, text="Export Accel Time History", command=quit)
    button_ex.config(height=2, width=23, state='disabled')
    button_ex.grid(row=crow, column=2, columnspan=2, padx=10, pady=3)

    button_ex_psd = tk.Button(top, text="Export PSD", command=quit)
    button_ex_psd.config(height=2, width=14, state='disabled')
    button_ex_psd.grid(row=crow, column=4, columnspan=1, padx=10, pady=3)

    ###############################################################################

    # start event-loop
    root.mainloop()

    return


def Natural_Clicked():
    if tk.BooleanVar.get(natural):
        NaturalEntry.grid(row=5, column=1)
        print('Checkbox 1 selected')
    else:
        NaturalEntry.grid_forget()
        print('Checkbox 1 unselected')


def Create_gif(filenames, gifname):
    frames = [imageio.imread(filename) for filename in filenames]
    imageio.mimsave(gifname, frames, format='GIF', duration=1)

    # Remove files
    for filename in set(filenames):
        os.remove(filename)


def PostProcessLoS(LoS_Dct):
    fr = []  # [1,1,1,1,2,2,2,2,3,3,3,3...]
    dx = []
    dy = []

    fr_2 = []  # [1,2,3,4...]
    avgX_lst = []
    avgY_lst = []

    for k in LoS_Dct.keys():
        ii = 0
        avgX = 0
        avgY = 0
        stdv = 0
        for v in LoS_Dct[k]:
            fr.append(k)
            dx.append(float((v[1])))
            dy.append(float((v[2])))
            ii += 1
            avgX += float((v[1]))
            avgY += float((v[2]))

        avgX = avgX / ii
        avgY = avgY / ii
        avgX_lst.append(avgX)
        avgY_lst.append(avgY)
        fr_2.append(k)

    # saving the results in path
    path = os.getcwd()
    print(path)
    LoS_dir_path = path + '\LoS_Results '
    # changing the cur working dir
    try:
        os.chdir(LoS_dir_path)
    except FileNotFoundError as e:
        print('Error! create a folder called "LoS_Results" manually or delete the OutSim folder and re-run'.format(e))

    # Plot:
    fig, axs = plt.subplots(nrows=2, figsize=(25, 10))
    fig.suptitle('LoS Results per Frame / Sub-Frame. (for urad: multiply by {})'.format(
        ParamsDict["micron_to_pixel"] / (ParamsDict["EFL"])))
    axs[0].plot(fr, dx, 'bo', label="LoS during integration time (Sub-Frame)")
    axs[0].plot(fr_2, avgX_lst, 'ro-', label="Avg. LoS per frame")
    axs[0].legend(loc='lower right')
    axs[0].set_xlabel("Frame #")
    axs[0].tick_params(axis='x', labelsize=5.5)
    axs[0].set_ylabel("DX [Pixel]")

    axs[1].plot(fr, dy, 'bo', label="LoS during integration time (Sub-Frame)")
    axs[1].plot(fr_2, avgY_lst, 'ro-', label="Avg. LoS per frame")
    axs[1].legend(loc='lower right')
    axs[1].set_xlabel("Frame #")
    axs[1].tick_params(axis='x', labelsize=5.5)
    axs[1].set_ylabel("DY [Pixel]")

    plt.savefig('LoS Results per Frame')

    # save Data
    import pandas as pd
    df = pd.DataFrame()
    IFOV = ParamsDict["micron_to_pixel"] / ParamsDict["EFL"]
    df["Frame"] = fr_2
    df["Avg_LoS_DX (pixel)"] = avgX_lst
    df["Avg_LoS_DY (pixel)"] = avgY_lst
    df["Avg_LoS_DX (urad)"] = df["Avg_LoS_DX (pixel)"] * IFOV
    df["Avg_LoS_DY (urad)"] = df["Avg_LoS_DY (pixel)"] * IFOV
    df.to_csv("Avg_LoS_Data.csv")
    #####################################################

    plt.figure("Hist LoS DX", figsize=(10, 10))
    # print("dx", dx, type(dx))
    dx = np.asarray(dx)
    dx = dx * (ParamsDict["micron_to_pixel"] / (ParamsDict["EFL"]))  # SWITCH TO URAD
    maxVal = dx.max()
    minVal = dx.min()
    mu0, std0 = norm.fit(dx)
    # print("maxVal", maxVal, minVal)
    plt.hist(dx, bins=51, density=True)
    xmin, xmax = plt.xlim()
    xx = np.linspace(xmin, xmax, 100)
    pp = norm.pdf(xx, mu0, std0)
    plt.plot(xx, pp, 'k', linewidth=2)

    plt.xlim(xmin=-maxVal, xmax=maxVal)
    plt.ylim(ymin=0, ymax=0.5)
    # plt.yscale('log')
    title0 = " Tilts Histogram of Simulated Movie; Fit Values: {:.2f} and {:.2f}".format(mu0, std0)
    plt.title(title0)
    plt.xlabel(' Angel (urad)')
    plt.ylabel(' Counts ')

    plt.savefig('Los_Pixels_Hist_X')
    ##############################################################

    plt.figure("Hist LoS DY", figsize=(10, 10))
    # print("dx", dx, type(dx))
    dy = np.asarray(dy)
    dy = dy * (ParamsDict["micron_to_pixel"] / (ParamsDict["EFL"]))  # SWITCH TO URAD
    maxVal = dy.max()
    minVal = dy.min()
    mu0, std0 = norm.fit(dy)
    # print("maxVal", maxVal, minVal)
    plt.hist(dy, bins=51, density=True)
    xmin, xmax = plt.xlim()
    xx = np.linspace(xmin, xmax, 100)
    pp = norm.pdf(xx, mu0, std0)
    plt.plot(xx, pp, 'k', linewidth=2)

    plt.xlim(xmin=-maxVal, xmax=maxVal)
    plt.ylim(ymin=0, ymax=0.5)
    # plt.yscale('log')
    title0 = " Tilts Histogram of Simulated Movie; Fit Values: {:.2f} and {:.2f}".format(mu0, std0)
    plt.title(title0)
    plt.xlabel(' Angel (urad)')
    plt.ylabel(' Counts ')

    plt.savefig('Los_Pixels_Hist_Y')

    # plt.show()

    return avgX_lst, avgY_lst


def get_fig_from_gauss_dist(dist_DX, dist_DY, samp_from_DX, samp_from_DY):
    # making probabilities vector probabilities = [0,1] based on DX,DY
    x_dx = np.linspace(dist_DX.ppf(.00001), dist_DX.ppf(.99999))
    x_dy = np.linspace(dist_DY.ppf(.00001), dist_DY.ppf(.99999))

    fig, axs = plt.subplots(2, 2, figsize=(10, 8))

    axs[0, 0].plot(x_dy, dist_DY.pdf(x_dy), color='k')
    axs[0, 0].axvline(samp_from_DY, color='r', ls='--')
    new_line = '\n'
    stats_str = f" z ={round(samp_from_DY, 3)} {new_line} PDF(z) ={round(dist_DY.pdf(samp_from_DY), 3)} {new_line} CDF(z) ={round(dist_DY.cdf(samp_from_DY), 3)}"
    x_lim_00 = axs[0, 0].get_xlim()
    y_lim = axs[0, 0].get_ylim()
    axs[0, 0].text(0.85 * x_lim_00[0], 0.8 * y_lim[1], stats_str, bbox=dict(facecolor='blue', alpha=0.1))
    axs[0, 0].set(
        title='PDF - DY'
    )

    axs[0, 1].plot(x_dy, dist_DY.cdf(x_dy))
    axs[0, 1].axhline(dist_DY.cdf(samp_from_DY), color='r', ls='--')
    axs[0, 1].set(
        title='CDF - DY'
    )

    ############ DX ##########################  :
    axs[1, 0].plot(x_dx, dist_DX.pdf(x_dx), color='k')
    axs[1, 0].axvline(samp_from_DX, color='r', ls='--')
    new_line = '\n'
    stats_str = f" z ={round(samp_from_DX, 3)} {new_line} PDF(z) ={round(dist_DX.pdf(samp_from_DX), 3)} {new_line} CDF(z) ={round(dist_DX.cdf(samp_from_DX), 3)}"
    x_lim_10 = axs[1, 0].get_xlim()
    y_lim = axs[1, 0].get_ylim()

    ## Set X limits to be euqal for the plots:
    max_x_lim = max(x_lim_00[1], x_lim_10[1])
    axs[0, 0].set_xlim(-1 * max_x_lim, max_x_lim)
    axs[0, 1].set_xlim(-1 * max_x_lim, max_x_lim)
    axs[1, 0].set_xlim(-1 * max_x_lim, max_x_lim)
    axs[1, 1].set_xlim(-1 * max_x_lim, max_x_lim)

    # retake
    x_lim_10 = axs[1, 0].get_xlim()
    y_lim = axs[1, 0].get_ylim()

    axs[1, 0].text(0.85 * x_lim_10[0], 0.8 * y_lim[1], stats_str, bbox=dict(facecolor='blue', alpha=0.1))
    axs[1, 0].set(
        title='PDF - DX'
    )

    axs[1, 1].plot(x_dx, dist_DX.cdf(x_dx))
    axs[1, 1].axhline(dist_DY.cdf(samp_from_DX), color='r', ls='--')

    axs[1, 1].set(
        title='CDF - DX'
    )
    return fig


def Plot_sample_from_gauss_dist(dist_DX, dist_DY, avgX_lst, avgY_lst):
    plt.close("all")
    filenames = []
    i = 1
    for samp_from_DX, samp_from_DY in zip(avgX_lst, avgY_lst):
        fig = get_fig_from_gauss_dist(dist_DX, dist_DY, samp_from_DX, samp_from_DY)
        filename = f'{i}.png'
        filenames.append(filename)
        plt.savefig(filename)
        plt.close()
        i = i + 1
        if i > 50:
            break

    ####### build gif #########################################
    gifname = "Show_Sampling.gif"
    Create_gif(filenames, gifname)
    return


def create_Out_folder():
    # Unpack necessary parameters:
    IntTime = ParamsDict["IntTime"]
    IntTime = IntTime / 1000  # ms to s
    cam_sample_rate = ParamsDict["cam_sample_rate"]
    H = ParamsDict["DetectorH"]
    V = ParamsDict["DetectorV"]
    ExagFactor = ParamsDict["ExagFactor"]

    image_folder_name = '\SimOutput '
    if natural:
        FrameInputStr = FrameInput.split(".")[0]
    else:
        FrameInputStr = "SynFrame"
    image_folder_name = image_folder_name + FrameInputStr + ' Exp %1.2fms Motion x%d' % (IntTime * 1000, ExagFactor)
    path = os.getcwd()
    Newdir_path = path + image_folder_name

    try:
        os.mkdir(Newdir_path)
        if PSNR_cb:
            PSNR_dir_path = Newdir_path + '\PSNRImages '
            os.mkdir(PSNR_dir_path)
        if LoS_on_detector_cb:
            LoS_dir_path = Newdir_path + '\LoS_Results '
            os.mkdir(LoS_dir_path)

    except OSError:
        print("Creation of the directory %s failed" % Newdir_path)
    os.chdir(Newdir_path)
    print(os.getcwd())


def save_PSNR_results(nominal_frame, blurred_frame, frame_number):
    path = os.getcwd()
    PSNR_dir_path = path + '\PSNRImages '
    # changing the cur working dir
    try:
        os.chdir(PSNR_dir_path)
    except FileNotFoundError as e:
        print('Error! create a folder called "PSNRImages" manually or delete the OutSim folder and re-run'.format(e))

    # "perfect PSNR"
    if frame_number == 1:
        PSNR = cv2.PSNR(nominal_frame, nominal_frame)
        print('PERFECT PSNR: {:.2f}dB'.format(PSNR))
        frame_title = '_Nominal Frame PSNR {:.2f} dB.jpg'.format(PSNR)
        cv2.imwrite(frame_title, nominal_frame)

    PSNR = cv2.PSNR(nominal_frame, blurred_frame)

    print('PSNR: {:.2f}dB'.format(PSNR))

    # # PSNR direct calculation
    # mse = np.mean((nominal_frame.astype(np.float64) / 255 - blurred_frame.astype(np.float64) / 255) ** 2)  # mean squared error
    # print('Own implementation: ', 10 * np.log10(1. / mse))

    frame_title = 'Frame #{} PSNR {:.2f} dB.jpg'.format(str(frame_number), PSNR)
    cv2.imwrite(frame_title, blurred_frame)

    # changing back the cur working dir
    os.chdir(path)
    return


def create_USAF_frame():
    # unpack necessary params

    H = ParamsDict["DetectorH"]
    V = ParamsDict["DetectorV"]
    pix_size = ParamsDict["micron_to_pixel"]  # micron
    Fnyq = round(1000 / (2 * pix_size), 3)
    MP = ParamsDict["EFL"] / ParamsDict["CFL"]  # Magnification Product
    MP = round(MP, 3)
    Exag = ParamsDict["ExagFactor"]
    IFOV = round(pix_size / ParamsDict["EFL"], 1)  # urad (um/m)

    ##Write parameters to clip scene
    txttowrite = "---SIMULATION PARAMETERS:--------------------------------------------\n"
    txttowrite += "Effective Focal Length: " + str(
        ParamsDict["EFL"] * 1000) + " [mm] " + "  |||   Collimator Focal Length: " + str(
        ParamsDict["CFL"] * 1000) + " [mm]" + "   |||   MP: x" + str(MP) + "\n"
    txttowrite += "Detector Size: " + str(ParamsDict["DetectorH"]) + " x " + str(
        ParamsDict["DetectorV"]) + " [pixels]" + "    |||   Pixel Size: " + str(
        ParamsDict["micron_to_pixel"]) + " [um]" + "   |||    IFOV.: " + str(IFOV) + " [urad]\n"
    txttowrite += "Frame rate: " + str(
        ParamsDict["cam_sample_rate"]) + " [FPS]" + "                 |||   Integration Time: " + str(
        ParamsDict["IntTime"]) + " [ms]" + "   |||    Nyquist Freq.: " + str(Fnyq) + " [cy/mm]\n"
    txttowrite += "(*Synthetic signal sampling rate: " + str(ParamsDict["SampRate"]) + " [Hz])\n"
    txttowrite += "----------------------------"

    # print(txttowrite)

    # Draws USAF Frame
    DataFrame = np.zeros((149, H, 3), np.uint8)
    new_frame = np.zeros((V, H, 3), np.uint8)
    DataFrame = cv2.copyMakeBorder(DataFrame, 0, 1, 0, 0, cv2.BORDER_ISOLATED, value=[255, 255, 255])

    ### Writing txt to DataFrame
    white = (255, 255, 255)
    # Write Params txt in Frame
    y0, dy = 15, 23
    for i, line in enumerate(txttowrite.split('\n')):
        y = y0 + i * dy
        cv2.putText(DataFrame, line, (10, y), cv2.FONT_HERSHEY_SIMPLEX, 0.4, white, 1)

    titletxt = "Resolution Target | Labeled in cycles per mm"
    if Exag > 1:
        titletxt += " // motion amplified x" + str(Exag)

    cv2.putText(DataFrame, titletxt, (int(H / 2) - 220, 130), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 255, 0), 1)

    credittxt = "Creator: Yarden Zaki"
    cv2.putText(DataFrame, credittxt, (int(H) - 135, 15), cv2.FONT_HERSHEY_SIMPLEX, 0.4, (0, 0, 255), 1)

    new_frame = Draw_USAF_Bars(new_frame)

    # Making a test of creating natural frame:  #### NATURAL
    if natural:
        ImgStr = os.path.dirname(os.getcwd())
        ImgStr = ImgStr + "/" + FrameInput
        print("ImgStr", ImgStr)
        new_frame = cv2.imread(ImgStr)
        hor = new_frame.shape[1]
        ver = new_frame.shape[0]
        center = (int(ver / 2), int(hor / 2))
        # new_frame[center[0] - int(V/2):center[0] + int(V/2), center[1] - int(H/2):center[1] + int(H/2)] = 0-
        new_frame = new_frame[center[0] - int(V / 2):center[0] + int(V / 2),
                    center[1] - int(H / 2):center[1] + int(H / 2)]
        # print(new_frame.shape)
        # cv2.imshow("ressss",new_frame)
        # cv2.waitKey(4000)
        # exit()

    if FF:
        # Save Nomional Image of Frame
        fname = "NominalFrame_" + str(ParamsDict["DetectorH"]) + " x " + str(ParamsDict["DetectorV"]) + "_" + str(
            1000 * ParamsDict["EFL"]) + "mm.jpg"
        NominalFrameToSave = cv2.vconcat([DataFrame, new_frame])
        cv2.imwrite(fname, NominalFrameToSave)
        cv2.imwrite("AnalyzedFrame.jpg", new_frame)

    return new_frame, DataFrame


def Draw_USAF_Bars(new_frame):
    # Unpack Necessary params:
    H = ParamsDict["DetectorH"]
    V = ParamsDict["DetectorV"]
    micron_to_pixel = ParamsDict["micron_to_pixel"]
    EFL = ParamsDict["EFL"]
    CFL = ParamsDict["CFL"]
    ShrinkFactor = EFL / CFL

    white = (255, 255, 255)
    pixel_to_mm = (1 / micron_to_pixel) * 1000

    frame_with_bars = new_frame.copy()

    SpatialFreqToPlot = [0.66, 1, 2, 3, 5, 6]  # cy/mmm

    xcoord = MarginOfFrame + 20
    ycoord = MarginOfFrame + 80
    dxx = 30
    max_y_coord = 0
    cv2.putText(frame_with_bars, "Cy / mm", (2, 10), cv2.FONT_HERSHEY_SIMPLEX, 0.35, white, 1)
    for SF in SpatialFreqToPlot:
        ResTxt = str(round(SF, 2))
        HalfCycle = (1 / SF) * pixel_to_mm * ShrinkFactor * 0.5
        # print("HalfCycle [pixels] , Spatial Freq  ",HalfCycle, SF)
        HalfCycle = int(round(float("%.4g" % HalfCycle), 0))  # square width - HALF CYCLE
        if 1 <= HalfCycle and xcoord + 12 * HalfCycle < H - MarginOfFrame - 20:
            # Drawing Rectangles
            for i in range(0, 6, 2):
                # print(i)
                # H - Bars
                frame_with_bars[ycoord:ycoord + 5 * HalfCycle,
                xcoord + i * HalfCycle:xcoord + HalfCycle * (i + 1)] = white
                # V - Bars
                frame_with_bars[ycoord + i * HalfCycle:ycoord + HalfCycle * (i + 1),
                xcoord + 6 * HalfCycle:5 * HalfCycle + xcoord + 6 * HalfCycle] = white
                # max y (so that the next line will be lower)
                max_y_coord = max(max_y_coord, ycoord + HalfCycle * (i + 1))
            # #DEBUG
            # #draw green dots over spatial samples
            # frame_with_bars[ycoord -3, xcoord] = (0, 255, 0)
            # frame_with_bars[ycoord -3, xcoord +HalfCycle] = (0, 255, 0)
            # frame_with_bars[ycoord -3, xcoord +2*HalfCycle] = (0, 255, 0)
            # print("position of pixels along spatial length (FULL CYCLE):",xcoord, xcoord +2*HalfCycle)
            # cv2.imshow("fbr1", frame_with_bars)

            # Drawing Target Txt:
            ResTxtSize = 0.9 * HalfCycle / 7
            if ResTxtSize < 0.4:
                ResTxtSize = 0.4
            if ResTxtSize > 2.2:
                ResTxtSize = 2.2
            # print("ResTxtSize",ResTxtSize)
            cv2.putText(frame_with_bars, ResTxt, (int(xcoord + 6 * HalfCycle - ResTxtSize * 50), ycoord - 20),
                        cv2.FONT_HERSHEY_SIMPLEX, ResTxtSize, white, 2)
            # Promoting index:
            xcoord = xcoord + 12 * HalfCycle + dxx

    NyqLimitFreq = (1000 / (2 * micron_to_pixel)) * ShrinkFactor
    JhonsonDetect = 2 * (pixel_to_mm * ShrinkFactor * 0.5) / 1.5  # HALF CYCLE!!
    JhonsonRegognize = 2 * (pixel_to_mm * ShrinkFactor * 0.5) / 6  # HALF CYCLE!!
    JhonsonIdentify = 2 * (pixel_to_mm * ShrinkFactor * 0.5) / 12  # HALF CYCLE!!
    SpatialFreqToPlot = [JhonsonIdentify, JhonsonRegognize, JhonsonDetect, NyqLimitFreq]  # cy/mm

    SF_Lst_Str = ["JhonsonIdentify ", "JhonsonRegognize", "JhonsonDetect   ", "NyqLimitFreq    "]
    kk = 0

    xcoord = MarginOfFrame + 20
    # print("ccccc",V,max_y_coord)
    if max_y_coord + 130 < V - MarginOfFrame - 20:
        ycoord = max_y_coord + 130
    else:
        ycoord = V - 50

    dxx = 30
    for SF in SpatialFreqToPlot:
        ResTxt = str(round(SF, 2))
        HalfCycle = (1 / SF) * pixel_to_mm * ShrinkFactor * 0.5
        HalfCycle = int(round(float("%.4g" % HalfCycle), 0))  # square width - HALF CYCLE
        print(SF_Lst_Str[kk], "HalfCycle [pixels] -----", HalfCycle, "FullCycle [pixels] -----", 2 * HalfCycle,
              "Spatial Freq ----  ", ResTxt)
        kk += 1
        if 1 <= HalfCycle and xcoord + 12 * HalfCycle < H - MarginOfFrame - 20:
            # Drawing Rectangles
            for i in range(0, 6, 2):
                # print(i)
                # H - Bars
                frame_with_bars[ycoord:ycoord + 5 * HalfCycle,
                xcoord + i * HalfCycle:xcoord + HalfCycle * (i + 1)] = white
                # V - Bars
                frame_with_bars[ycoord + i * HalfCycle:ycoord + HalfCycle * (i + 1),
                xcoord + 6 * HalfCycle:5 * HalfCycle + xcoord + 6 * HalfCycle] = white

                # #DEBUG
                # #draw green dots over spatial samples
                # frame_with_bars[ycoord -3, xcoord] = (0, 255, 0)
                # frame_with_bars[ycoord -3, xcoord +HalfCycle] = (0, 255, 0)
                # frame_with_bars[ycoord -3, xcoord +2*HalfCycle] = (0, 255, 0)
                # print("position of pixels along spatial length (FULL CYCLE):",xcoord, xcoord +2*HalfCycle)
                # cv2.imshow("fbr1", frame_with_bars)

                # Drawing Target Txt:
            ResTxtSize = 0.9 * HalfCycle / 7
            if ResTxtSize < 0.4:
                ResTxtSize = 0.4
            if ResTxtSize > 2.2:
                ResTxtSize = 2.2
            # print("ResTxtSize", ResTxtSize)
            cv2.putText(frame_with_bars, ResTxt, (int(xcoord + 6 * HalfCycle - ResTxtSize * 50), ycoord - 20),
                        cv2.FONT_HERSHEY_SIMPLEX, ResTxtSize, white, 2)
            # Promoting index:
            xcoord = xcoord + 12 * HalfCycle + dxx

    return frame_with_bars


def Check_corr_between_two_lists(l1, l2):
    a1 = l1
    a2 = l2

    print(np.corrcoef(a1, a2))

    return


def BlurBasedOnIntTime(new_frame, n_exposure, data_dx, ECDF_data_dx, sorted_data_dx, nn, data_dy, sorted_data_dy,
                       ECDF_data_dy, frame_num):
    global FF
    x, y, _ = new_frame.shape
    canvas = new_frame.copy()

    # #tmppp
    # cv2.imshow("canvas0",canvas)
    # canvas[MarginOfFrame:x-MarginOfFrame,MarginOfFrame:y-MarginOfFrame]=(130,130,130)
    # cv2.imshow("canvas1", canvas)
    # cv2.waitKey(5000)

    if FF:
        # video Writer for making a movie of what is happening during the int time -------------------------
        fourcc = cv2.VideoWriter_fourcc(*'MJPG')
        fps_vid_writer = int(7)

        # Unpack necessary parameters:
        IntTime = ParamsDict["IntTime"]
        IntTime = IntTime / 1000  # ms to s
        ExagFactor = ParamsDict["ExagFactor"]
        SampRate = ParamsDict["SampRate"]
        vidFname = "Sampling "
        vidFname = vidFname + '%dHz signal during exposure %1.2fms Motion x%d.avi' % (
        SampRate, IntTime * 1000, ExagFactor)
        vidWriter = cv2.VideoWriter(vidFname, fourcc, fps_vid_writer, (int(y), int(x)), True)

    LoS_results = []
    ii = 1
    """
    Main physics loop - Features Added:
    1. Matching the CDF of DX and DY to promise correlation between them on high amp. vibrations
    2.  during "n_exposure" (exposure time) - if the first random sampling is with CDF "Phi" - the next frames will be based on its neighbors in TH
        (find argwhere of data_dx==samp_from_DX) then arg = arg +1 for consecutive "exposures". DY works with the same logic

    """
    for fr in range(0, n_exposure):  # The detector samples all of the "Analog" Data during the exposure time....
        # vibrating the area around the object: --------------------------------------
        # Calc the delta movements -- based on CDF - see C:\Users\E010236\Desktop\Python\LoS_Simulation\roy ayash\plot_dist_python
        if fr == 0:
            samp_from_DX = np.random.choice(data_dx)  # randomly sampled datapoint from data_dx dist.
            indx_sampld_DX = np.argwhere(data_dx == samp_from_DX)[0][0]
        else:
            try:
                indx_sampld_DX = indx_sampld_DX + 1
                samp_from_DX = data_dx[indx_sampld_DX]
            except IndexError:
                indx_sampld_DX = indx_sampld_DX - 1
                samp_from_DX = data_dx[indx_sampld_DX]

        CDF = round(ECDF_data_dx(samp_from_DX), 3)  # calculating the CDF of this data point
        # sampling from SORT (SORTED DY) based on the calculated CDF
        if fr == 0:
            try:
                samp_from_DY = sorted_data_dy[int(nn * CDF)]
                indx_sampld_DY = np.argwhere(data_dy == samp_from_DY)[0][0]
            except IndexError:
                if CDF > 0.5:
                    CDF = CDF - 0.01
                else:
                    CDF = CDF + 0.01
                samp_from_DY = sorted_data_dy[int(nn * CDF)]
                indx_sampld_DY = np.argwhere(data_dy == samp_from_DY)[0][0]

        else:
            try:
                indx_sampld_DY = indx_sampld_DY + 1
                samp_from_DY = data_dy[indx_sampld_DY]
            except IndexError:
                indx_sampld_DY = indx_sampld_DY - 1
                samp_from_DY = data_dy[indx_sampld_DY]

        # # for validation only:  - MUST TO UNCOMMENT
        # c1 = CDF
        # ECDF_data_dy = ECDF(sorted_data_dy)
        # c2 = round(ECDF_data_dy(samp_from_DY), 3)
        # print("EQUAL CDFS",c1==c2)
        # ####---------------------------------------

        dx = int(round(float("%.4g" % samp_from_DX), 0))  # Currently ROUNDING the dx,dy - think about improvements
        dy = int(round(float("%.4g" % samp_from_DY), 0))
        if LoS_on_detector_cb:
            LoS_results.append([ii, "%.4g" % samp_from_DX, "%.4g" % samp_from_DY])

        # print("counting...", ii , " / " , n_exposure)
        ii = ii + 1
        # print("dx","dy, index                           ",dx, dy ,startindx+fr)

        # Crop the ROI based on the wanted BoxSize
        cropped_area = new_frame[MarginOfFrame:x - MarginOfFrame, MarginOfFrame:y - MarginOfFrame]

        # creating result frame
        tmp_frame = canvas.copy()
        # FF moving is the "new stationary"
        moving_crop = cropped_area.copy()
        # pasting the or result to tmp_frame
        try:
            tmp_frame[MarginOfFrame + dy:x - MarginOfFrame + dy,
            MarginOfFrame + dx:y - MarginOfFrame + dx] = moving_crop.copy()
        except ValueError as e:
            print('Error! Frame margin of {} is too small. Set a greater value'.format(MarginOfFrame))

        # cv2.imshow("tmp_frame", tmp_frame)

        if fr == 0:
            new_line = "\n"
            to_print = f" sampled DX , DY : {dx} {dy} {new_line} CDF of sampled  DX ,DY: {CDF} {new_line}"
            print(to_print)

            # initializing out frame
            BlurredOutFrame = tmp_frame.copy()

        else:  # if fr >0
            alpha = 1.0 / (ii - 1)
            beta = 1.0 - alpha
            BlurredOutFrame = cv2.addWeighted(tmp_frame, alpha, BlurredOutFrame, beta, 0.0)
            ### Black border padding
            BlurredOutFrame[0:, 0:MarginOfFrame] = 0
            BlurredOutFrame[0:, -1 * MarginOfFrame:] = 0
            BlurredOutFrame[0:MarginOfFrame, 0:] = 0
            BlurredOutFrame[-1 * MarginOfFrame:, 0:] = 0
            # BlurredOutFrame = cv2.copyMakeBorder(BlurredOutFrame,MarginOfFrame,MarginOfFrame,MarginOfFrame,
            #                                    MarginOfFrame,cv2.BORDER_CONSTANT,value=0)

            # print("ourframeSHAPE",BlurredOutFrame.shape)

        if debug:
            # This is what happened during the exposure time on the detector
            cv2.imshow("BlurredOutFrame_tmp", BlurredOutFrame)
            cv2.waitKey(20) & 0xFF
        if FF:
            vidWriter.write(BlurredOutFrame)

            cv2.waitKey(5)

    # Adding Los results to LoS_Dict.
    LoS_Dct[str(frame_num)] = LoS_results

    # reducing brightness level to get more realistic outcome
    # BrightFactor = 1
    # if natural:
    #     BrightFactor = 1
    # BlurredOutFrame [MarginOfFrame + dy:x - MarginOfFrame + dy,MarginOfFrame + dx:y - MarginOfFrame + dx]  = BrightFactor * BlurredOutFrame [MarginOfFrame + dy:x - MarginOfFrame + dy,MarginOfFrame + dx:y - MarginOfFrame + dx]

    # BGR
    BGRframe = BlurredOutFrame.copy()
    # BGR to GRAY
    BlurredOutFrame = cv2.cvtColor(BlurredOutFrame, cv2.COLOR_BGR2GRAY).copy()

    # Releasing the VidWriter
    if FF:
        vidWriter.release()

    # plotting the "working area"
    if debug:
        WRframe = BlurredOutFrame.copy()
        WRframe[MarginOfFrame + dy:x - MarginOfFrame + dy, MarginOfFrame + dx:y - MarginOfFrame + dx] = 255
        # cv2.imshow("WRframe",WRframe)

    FF = False
    return BlurredOutFrame, BGRframe


def create_blured_frame(sample_rate_of_input, dur, input_dx, input_dy):
    # Unpack necessary parameters:
    IntTime = ParamsDict["IntTime"]
    IntTime = IntTime / 1000  # ms to s
    cam_sample_rate = ParamsDict["cam_sample_rate"]
    H = ParamsDict["DetectorH"]
    V = ParamsDict["DetectorV"]
    ExagFactor = ParamsDict["ExagFactor"]

    # add sampling rate to ParamsDict
    ParamsDict["SampRate"] = round(sample_rate_of_input, 1)

    # Total time:
    SamplingIncrement = int(
        sample_rate_of_input / cam_sample_rate)  # 3000Hz signal samples with 30fps camera and therefore the signal should be sampled at increments of 100 (the factor)

    SamplingIncrementDuringExposure = int(
        sample_rate_of_input / (
                    1 / IntTime))  # 3000Hz signal samples with 1/int_time detector and therefore the signal should be sampled at increments of sample_rate_of_input*IntTime (the factor)

    duration = dur  # [s]
    # Total N
    n = cam_sample_rate * duration
    n_exposure = SamplingIncrementDuringExposure  # number of samples during the exposure time (integration time)

    # video Writer -------------------------
    fourcc = cv2.VideoWriter_fourcc(*'MJPG')
    fps_vid_writer = int(cam_sample_rate)
    vidFname = "USAF "
    vidFname = vidFname + '%dsec %dFPS Exposure %1.2fms Motion x%d.avi' % (
    duration, fps_vid_writer, IntTime * 1000, ExagFactor)
    vidWriter = cv2.VideoWriter(vidFname, fourcc, fps_vid_writer, (H, V + 150), True)
    # ---- BGR video writer
    vidFname2 = "BGR " + 'exposure %1.2fms Motion x%d.avi' % (IntTime * 1000, ExagFactor)
    vidWriter2 = cv2.VideoWriter(vidFname2, fourcc, fps_vid_writer, (H, V), True)
    # -----------------------------------------

    data_dx = input_dx * ExagFactor
    data_dy = input_dy * ExagFactor

    # Empirical CDF of DX - our real array is not pure gaussian. - CONSIDER DO IT ONCE OUT OF LOOP
    ECDF_data_dx = ECDF(data_dx)
    ECDF_data_dy = ECDF(data_dy)
    nn = len(data_dy)
    sorted_data_dy = sorted(data_dy)
    sorted_data_dx = sorted(data_dx)

    ii = 1
    new_frame, DataFrame = create_USAF_frame()

    for fr in range(0, n):
        print("counting...", ii, " / ", n)
        k = ii
        LoS_Dct[str(k)] = None
        ii = ii + 1
        # vibrating the area around the object: --------------------------------------
        x, y, _ = new_frame.shape

        if True:
            print("Number of Sub-Frames during the int Time: ", n_exposure)
            # Blurring Based on exposure time (physical method)

            BlurResults = BlurBasedOnIntTime(new_frame, n_exposure, data_dx, ECDF_data_dx, sorted_data_dx, nn, data_dy,
                                             sorted_data_dy, ECDF_data_dy, k)
            OutFrame, BGRoutframe = BlurResults[0], BlurResults[1]
            OutFrame = cv2.cvtColor(OutFrame, cv2.COLOR_GRAY2BGR)
            if PSNR_cb:
                save_PSNR_results(new_frame, BGRoutframe, k)
            vidWriter2.write(BGRoutframe)

        # print(type(OutFrame), type(DataFrame))
        # print(OutFrame.shape, DataFrame.shape)
        OutFrame = cv2.vconcat([DataFrame, OutFrame])
        # OutFrame = cv2.cvtColor(OutFrame.copy(), cv2.COLOR_GRAY2BGR).copy()
        cv2.imshow('OutFrame', OutFrame)
        cv2.waitKey(50) & 0xFF
        vidWriter.write(OutFrame)
        # print("fr",fr)

    cv2.waitKey(600) & 0xFF
    cv2.destroyAllWindows()
    vidWriter.release()
    vidWriter2.release()
    return


def SimulateLoS():
    global EFL, CFL, DetectorH, DetectorV, micron_to_pixel, cam_sample_rate, IntTime, ExagFactor
    global debug, MarginOfFrame, natural, FrameInput
    global PSNR_cb, LoS_on_detector_cb
    global FF
    global LoS_Dct, ParamsDict
    LoS_Dct = {}  # With a structure of: key: frame_number (str) value: [[sub-frame num during int time (int) , los_dx_pixels (str), los_dy_pixels (str) ],,... ]
    FF = True

    # Collecting Data from user input: ----------------
    FrameInput = tk.StringVar.get(FrameInput)
    MarginOfFrame = tk.IntVar.get(MarginOfFrame)
    debug = tk.BooleanVar.get(debug)
    PSNR_cb = tk.BooleanVar.get(PSNR_cb)
    LoS_on_detector_cb = tk.BooleanVar.get(LoS_on_detector_cb)
    natural = tk.BooleanVar.get(natural)
    # ---
    ParamsDict = {}
    ParamsDict["EFL"] = float(tk.StringVar.get(EFL)) / 1000  # [m]
    ParamsDict["CFL"] = float(tk.StringVar.get(CFL)) / 1000  # [m]
    ParamsDict["DetectorH"] = int(float(tk.StringVar.get(DetectorH)))
    ParamsDict["DetectorV"] = int(float(tk.StringVar.get(DetectorV)))
    ParamsDict["micron_to_pixel"] = float(tk.StringVar.get(micron_to_pixel))
    ParamsDict["cam_sample_rate"] = int(float(tk.StringVar.get(cam_sample_rate)))
    ParamsDict["IntTime"] = float(tk.StringVar.get(IntTime))
    ParamsDict["ExagFactor"] = float(tk.StringVar.get(ExagFactor))
    for v in ParamsDict.values():
        print(v, type(v))
    # ------------------------------------------------------------------

    # Create output folder
    create_Out_folder()

    if not debug:
        tkMessageBox.showinfo("Info", "Enter Elevation LoS error PSD - 2 columns txt file [Hz, rad^2/Hz]")
        T_el, TH_el = psd_syn.psd_syn_main()
        TH_dy_detector_pix = TH_el * ParamsDict["EFL"] * 1e6 / ParamsDict["micron_to_pixel"]  # See below
        # TH [rad]  * EFL [m] = dy[m] -> multiply by 1e6 -> dy [microns] -> divide by micron_to_pixel -> dy[pixels]

        tkMessageBox.showinfo("Info", "Enter Azimuth LoS error PSD - 2 columns txt file [Hz, rad^2/Hz]")
        T_az, TH_az = psd_syn.psd_syn_main()
        TH_dx_detector_pix = TH_az * ParamsDict["EFL"] * 1e6 / ParamsDict["micron_to_pixel"]

        ###################################
        print(" ")
        print("Signal stats DY [pixels]:")
        sr_DY, dt, mean, sd, rms, skew, kurtosis, dur0 = signal_stats(T_el, TH_dy_detector_pix, len(T_el))

        print(" ")
        print("Signal stats DX [pixels]:")
        sr_DX, dt1, mean1, sd1, rms1, skew1, kurtosis1, dur1 = signal_stats(T_az, TH_dx_detector_pix, len(T_az))

        dur = int(min(dur0, dur1) - 1)
        ###################
        print(" ")
        print("Writing DX , DY time history .... ")
        a_output_file_path = "Out Time History " + "_DY_Pixels" + ".txt"
        file_path = a_output_file_path.rstrip('\n')
        nn = len(TH_dy_detector_pix) - 2
        print("nn", nn)
        WriteData2(nn, T_el, TH_dy_detector_pix, file_path)
        # ---
        a_output_file_path = "Out Time History " + "_DX_Pixels" + ".txt"
        file_path = a_output_file_path.rstrip('\n')
        nn = len(TH_dx_detector_pix) - 2
        WriteData2(nn, T_az, TH_dx_detector_pix, file_path)

    if debug:
        ################TEMPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
        filepathdx = 'C:/Users/E010236/Desktop/Python/LoS_Simulation/debug_folder/Out Time History _DX_Pixels.txt'
        filepathdy = 'C:/Users/E010236/Desktop/Python/LoS_Simulation/debug_folder/Out Time History _DY_Pixels.txt'
        import os
        FileNamedx = os.path.split(filepathdx)

        print("FileNamedx", FileNamedx)
        FileNamedx = FileNamedx[-1]
        FileNamedx = FileNamedx.split('.')
        FileNamedx = FileNamedx[0]
        print("FileNamedx", FileNamedx)

        #
        if not os.path.exists(filepathdx):
            print("This file doesn't exist")
        #
        if os.path.exists(filepathdx):
            print("This file exists")
            print(" ")
            infile = open(filepathdx, "rb")
            lines = infile.readlines()
            infile.close()

            a = []
            b = []
            num = 0
            for line in lines:
                #
                if sys.version_info[0] == 3:
                    line = line.decode(encoding='UTF-8')

                if re.search(r"(\d+)", line):  # matches a digit
                    iflag = 0
                else:
                    iflag = 1  # did not find digit
                #
                if re.search(r"#", line):
                    iflag = 1
                #
                if iflag == 0:
                    line = line.lower()
                    if re.search(r"([a-d])([f-z])", line):  # ignore header lines
                        iflag = 1
                    else:
                        line = line.replace(",", " ")
                        col1, col2 = line.split()
                        a.append(float(col1))
                        b.append(float(col2))
                        num = num + 1

            T_az = np.array(a)
            TH_dx_detector_pix = np.array(b)

            FileNamedy = os.path.split(filepathdy)

            print("FileNamedy", FileNamedy)
            FileNamedy = FileNamedy[-1]
            FileNamedy = FileNamedy.split('.')
            FileNamedy = FileNamedy[0]
            print("FileNamedy", FileNamedy)

            #
        if not os.path.exists(filepathdy):
            print("This file doesn't exist")
        #
        if os.path.exists(filepathdy):
            print("This file exists")
            print(" ")
            infile = open(filepathdy, "rb")
            lines = infile.readlines()
            infile.close()

            a = []
            b = []
            num = 0
            for line in lines:
                #
                if sys.version_info[0] == 3:
                    line = line.decode(encoding='UTF-8')

                if re.search(r"(\d+)", line):  # matches a digit
                    iflag = 0
                else:
                    iflag = 1  # did not find digit
                #
                if re.search(r"#", line):
                    iflag = 1
                #
                if iflag == 0:
                    line = line.lower()
                    if re.search(r"([a-d])([f-z])", line):  # ignore header lines
                        iflag = 1
                    else:
                        line = line.replace(",", " ")
                        col1, col2 = line.split()
                        a.append(float(col1))
                        b.append(float(col2))
                        num = num + 1

            T_el = np.array(a)
            TH_dy_detector_pix = np.array(b)

        print(" ")
        print("Signal stats DY [pixels]:")
        sr_DY, dt, mean, sd, rms, skew, kurtosis, dur0 = signal_stats(T_el, TH_dy_detector_pix, len(T_el))

        print(" ")
        print("Signal stats DX [pixels]:")
        sr_DX, dt1, mean1, sd1, rms1, skew1, kurtosis1, dur1 = signal_stats(T_az, TH_dx_detector_pix, len(T_az))

        dur = int(min(dur0, dur1) - 1)

        ############ END TEMPPPPPPPPPPPPPPPPP

    # Check_corr_between_two_lists(TH_dx_detector_pix, TH_dy_detector_pix)
    # exit()

    print(ParamsDict)

    # getting data_dx , data_dy stats
    mu1, sigma1 = np.mean(TH_dx_detector_pix), np.std(TH_dx_detector_pix)  # DX stats:
    mu2, sigma2 = np.mean(TH_dy_detector_pix), np.std(TH_dy_detector_pix)  # DY stats:

    # Making a perfect normal distribuitions with DX,DY stats:
    Gauss_dist_DX = sps.norm(loc=mu1, scale=sigma1)
    Gauss_dist_DY = sps.norm(loc=mu2, scale=sigma2)

    create_blured_frame(sr_DX, dur, TH_dx_detector_pix, TH_dy_detector_pix)
    if LoS_on_detector_cb:
        avgX_lst, avgY_lst = PostProcessLoS(LoS_Dct)
        Plot_sample_from_gauss_dist(Gauss_dist_DX, Gauss_dist_DY, avgX_lst, avgY_lst)


#########################################-------------------------------------------------------------------------

# create root window
global debug, MarginOfFrame, natural
root = tk.Tk()

w, h = root.winfo_screenwidth(), root.winfo_screenheight()
w = int(2. * (w * 0.18))  # do not change
h = int(2. * (h * 0.4))  # do not change
root.geometry("%dx%d+0+0" % (w, h))

root.title("LoS Jitter Simulation & Visualization ver 1.0  by Yarden Zaki")

###############################################################################

crow = 1

hwtext1 = tk.Label(root, text='Simulates the effect of LoS error on known EO systems')
hwtext1.grid(row=crow, column=0, columnspan=3, padx=8, pady=7, sticky=tk.W)

crow += 1

hwtext1 = tk.Label(root, text='Note: bla bla')
hwtext1.grid(row=crow, column=0, columnspan=3, padx=8, pady=5, sticky=tk.W)

###############################################################################

crow += 1

hwtext2 = tk.Label(root, text='#' * 20)
hwtext2.grid(row=crow, column=0, columnspan=1, padx=8, pady=8)

hwtext2 = tk.Label(root, text='-' * 20)
hwtext2.grid(row=crow, column=1, columnspan=1, padx=8, pady=8)

###############################################################################
crow += 1
debug = tk.BooleanVar(root, False)
tk.Checkbutton(root, text="Debug Mode", variable=debug).grid(row=crow, column=0, sticky="w")

crow += 1
natural = tk.BooleanVar(root, False)
tk.Checkbutton(root, text="Simulate Natural Frame", variable=natural, command=Natural_Clicked).grid(row=crow, column=0,
                                                                                                    sticky="w")
FrameInput = tk.StringVar(root, value='Input_Frame.png')
NaturalEntry = tk.Entry(root, width=30, textvariable=FrameInput)

crow += 1
MarginOfFrame = tk.IntVar(root, 15)
tk.Label(root, text='Frame Margin [pixels]').grid(row=crow, column=0, sticky="w")
tk.Entry(root, textvariable=MarginOfFrame, width=5).grid(row=crow, column=1, sticky="w")

###############################################################################
crow += 1
hwtext2 = tk.Label(root, text='Simulation Parameters:  ' + '-' * 10).grid(row=crow, column=0, columnspan=1, pady=8,
                                                                          sticky="w")
###############################################################################
crow += 1
EFL = tk.StringVar(root, "1000")
tk.Label(root, text='Effective Focal Length of the System [mm]').grid(row=crow, column=0, sticky="w")
tk.Entry(root, textvariable=EFL, width=5).grid(row=crow, column=1, sticky="w")

crow += 1
CFL = tk.StringVar(root, "3000")
tk.Label(root, text='Effective Focal Length of the Collimator [mm]').grid(row=crow, column=0, sticky="w")
tk.Entry(root, textvariable=CFL, width=5).grid(row=crow, column=1, sticky="w")

crow += 1
DetectorH = tk.StringVar(root, "1920")
tk.Label(root, text='Detector Size Horizontal [pixels] ').grid(row=crow, column=0, sticky="w")
tk.Entry(root, textvariable=DetectorH, width=5).grid(row=crow, column=1, sticky="w")

crow += 1
DetectorV = tk.StringVar(root, "1080")
tk.Label(root, text='Detector Size Vertical [pixels]').grid(row=crow, column=0, sticky="w")
tk.Entry(root, textvariable=DetectorV, width=5).grid(row=crow, column=1, sticky="w")

crow += 1
micron_to_pixel = tk.StringVar(root, "5.86")
tk.Label(root, text='Pixel Size [um]').grid(row=crow, column=0, sticky="w")
tk.Entry(root, textvariable=micron_to_pixel, width=5).grid(row=crow, column=1, sticky="w")

crow += 1
cam_sample_rate = tk.StringVar(root, "25")
tk.Label(root, text='Camera Sample Rate [fps]').grid(row=crow, column=0, sticky="w")
tk.Entry(root, textvariable=cam_sample_rate, width=5).grid(row=crow, column=1, sticky="w")

crow += 1
IntTime = tk.StringVar(root, "5")
tk.Label(root, text='Integration Time [ms]').grid(row=crow, column=0, sticky="w")
tk.Entry(root, textvariable=IntTime, width=5).grid(row=crow, column=1, sticky="w")

crow += 1
ExagFactor = tk.StringVar(root, "1")
tk.Label(root, text='Motion Magnification Factor').grid(row=crow, column=0, sticky="w")
tk.Entry(root, textvariable=ExagFactor, width=5).grid(row=crow, column=1, sticky="w")
#################################################################################
crow += 1
hwtext2 = tk.Label(root, text='Output Selection:  ' + '-' * 10).grid(row=crow, column=0, columnspan=1, pady=8,
                                                                     sticky="w")
crow += 1
PSNR_cb = tk.BooleanVar(root, False)
tk.Checkbutton(root, text="PSNR", variable=PSNR_cb).grid(row=crow, column=0, sticky="w")
crow += 1
LoS_on_detector_cb = tk.BooleanVar(root, False)
tk.Checkbutton(root, text="LoS Data [pix]", variable=LoS_on_detector_cb).grid(row=crow, column=0, sticky="w")

crow += 1
buttonLoad = tk.Button(root, text='Load PSD results', command=Call_load_psd)
buttonLoad.grid(row=crow, column=0, padx=4, pady=30)
buttonLoad.config(height=2, width=30)

buttonStart = tk.Button(root, text='Start Simulation', command=SimulateLoS)
buttonStart.grid(row=crow, column=1, padx=4, pady=30)
buttonStart.config(height=2, width=20)

button_quit = tk.Button(root, text="Quit", command=lambda root=root: quit(root))
button_quit.grid(row=crow, column=2, padx=4, pady=30)
button_quit.config(height=2, width=20)

###############################################################################

# start event-loop

root.mainloop()
