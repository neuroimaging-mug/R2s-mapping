# 2D gradient-echo modelling 

This repository cotains MATLAB scripts used in Soellrad et al. 2019 for R2* and myelin water fraction (MWF) estimation in presence of macroscopic field variations for 2D spoiled gradient-echo sequences. 

To demonstrate the usage one simulation example and two in-vivo examples of a single subject for R2* and MWF mapping are provided. The in-vivo examples contain all postprocessing steps from coil combination to the R2* and MWF maps obtained with the different models. 

If you find the provided scripts useful for your work, please cite: 

Soellradl M, Lesch A, Strasser J, et al. 
  Assessment and correction of macroscopic field variations in 2D spoiled gradient-echo sequences. Magn Reson Med. 2019;00:1–14.
  https ://doi.org/10.1002/mrm.28139

Please also cite the [numerical Bloch solver](https://github.com/IMTtugraz/rfcontrol) used for slice profile estimation:

Aigner CS, Clason C, Rund A, Stollberger R. 
  Efficient high-resolution RF pulse design applied to simultaneous multi-slice excitation.
  J Magn Reson. 2016;263:33–44.



## Dependencies 

The scripts were tested with [MATLB R2016b](https://www.matlab.com). Besides the external tools, which are provided in the repository, [FSL 5.0.9](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki) is required for PRELUDE phase unwrapping and brain extraction with BET to process data as described in the data. If FSL is not installed different options are available to perform calculations without FSL.   

## Examples

### Example 1: Forward simulation
The example shows the effect of the Glsice polarity on the dephasing of the signal in presence of a macroscopic field gradient Gz, which was used to generate Figure 3 in the paper. 

### Example 2: R2s mapping <a href="https://doi.org/10.5281/zenodo.3600319"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.3600319.svg" alt="DOI"></a>
In this example R2s estimated with the proposed models. Before running the example, the input data and the results (optional) have to be downloaded by clicking on the zenodo link. 

### Example 3: MWF mapping <a href="https://doi.org/10.5281/zenodo.3600319"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.3600319.svg" alt="DOI"></a>

MWF maps with the different models are estimated as described in the paper for a single subject. In case the input data has not be downloaded please download input and results (optional) from the zenodo link. 

## Authors

* **Martin Soellradl** - *Initial work* - [mSoell](https://github.com/mSoell)


## License

This software is published under GNU GPLv3. 



