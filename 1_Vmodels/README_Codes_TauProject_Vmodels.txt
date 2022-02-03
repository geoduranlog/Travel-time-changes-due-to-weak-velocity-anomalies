======================= CODES Tau Project - Vmodels =====================
A. Duran 
Feb 2021

Here you find the codes for the Tau Project regarding the velocity models that will be put as inputs in the wave propagation simulations using SPECFEM2D


1. random_model_elastic_LargeSizeProduction_sigma20.m
It creates several realisations of 2D numerical velocity models using the Von-Karman autocorrelation function in the format used for the SPECFEM2D software. You specify the parameters: correlation length (ax, ay), density, etc

In this case sigma=0.2 indicates that the model has velocity fluctuations of 20% with respect to the reference (vref -> vp).

The folder "Models_example"  contains an example of how the model's file looks like.


2. Modify_Vmodel_at_rprime.m
This code is useful if you want to place a perturbation in the models and tests how the waveforms change before and after the perturbation, or also to test several perturbation strengths.  The code simply put a velocity perturbation (i.e., dvp, dvs, or both at the same time) over a fixed point in space (called "rprime") and export the models in a SPECDEM2D format.   

