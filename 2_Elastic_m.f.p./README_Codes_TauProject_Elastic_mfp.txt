======================= CODES Tau Project - Elastic mean free paths =====================
A. Duran 
Feb 11 2021

Here you find the codes to characterise the numerical media, i.e., computing the transport and elastic mean free paths


The synthetic waveforms are obtained from simulations on SPECFEM2D over 40 realisations of numerical media. The source is placed in the center of the medium and several stations are placed following a cross-pattern with the source in the center.
The waveforms are saved on the file “wave_coh_Model1-40.mat”, the recorded waveform is  the displacement.

——To get Scatt Mean free Paths——
1. A1_SMFP_coherent: Get the linear regression. Use P-arrival time as your start for the time window tp. I used a time window of tp:tp+tw    ->  Then you’re focusing on the 1st arrival which is the coherent part of the waveform.



——To get Transport Mean free Paths——

2. B1_TMFP_incoherent_diffusionfit:  Get the slide window avg numerical Intensities and the Diffusion coefficient you want to test.  In our case, we don’t test attenuation values.

3. B2_FIT_visu_incoherent.m:  Test Diffussion coefficient on Itheo and compare with Inum avg. Then Take the average of the Diff. coefficient and get the transport m.f.p.







