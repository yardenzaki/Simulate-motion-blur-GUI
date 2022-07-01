# Simulate-motion-blur-GUI
 Creating a simulation video from frame and LoS error analysis results

![img.png](img.png)

The desired ROI will vibrate according to the LoS error and statistical analysis
The analysis neglect any high - order image artifacts (defocus,astigmatism,coma etc..) and simulates only rigid body tilts ( first order Zernike polynoms)
How does it work?
1. Getting tilts LoS PSD Data [rad^2/Hz vs. Hz].
2. Converting to displacements LoS PSD with w.r.t the detector plane by multiplying with the (focal len)^2
3. Converting it to time history equivalent tilts [rad vs. time].
4. Render the corresponding image motion based on the results.


Dependencies:
* import psd_syn
* from tompy import enter_float, signal_stats, WriteData2