# BSReconFramework
Framework for the reconstruction of B1+ maps out of highly subsampled Bloch Siegert data

## Contributors
* Andreas Lesch (Graz University of Technology)
* [Matthias Schloegl](http://www.tugraz.at/institute/imt/people/schloegl/) (Graz University of Technology)
* [Martin Holler](http://imsc.uni-graz.at/hollerm) (University of Graz) 
* [Kristian Bredies](http://imsc.uni-graz.at/bredies) (University of Graz) 

## License
This software is published under GNU GPLv3. In particular, all source code is provided "as is" without warranty of any kind, either expressed or implied. For details, see the attached LICENSE.

## General Information
This framework is an open source software for the reconstruciton of B1+ maps out of highly subsampled Bloch Siegert data. The current version offers a CPU implementation in MATLAB
and a GPU implementation based on CUDA and the reconstruciton framework [AVIONIC](https://github.com/IMTtugraz/AVIONIC.git). 

If you use this software please cite:
* Schloegl, Matthias and Holler, Martin and Schwarzl, Andreas and Bredies, Kristian and Stollberger, Rudolf. <br>
  __Infimal convolution of total generalized variation functionals for dynamic MRI.__<br>
  _Magnetic Resonance in Medicine_, 2017; 78(1):142-155<br>
  doi: [10.1002/mrm.26352](http://onlinelibrary.wiley.com/doi/10.1002/mrm.26352/full) <br>
  OA: [EuroPubMed](http://europepmc.org/articles/PMC5553112)

```
@article {MRM:MRM26352,
author = {Schloegl, Matthias and Holler, Martin and Schwarzl, Andreas and Bredies, Kristian and Stollberger, Rudolf},
title = {Infimal convolution of total generalized variation functionals for dynamic MRI},
journal = {Magnetic Resonance in Medicine},
volume = {78},
number = {1},
issn = {1522-2594},
url = {http://dx.doi.org/10.1002/mrm.26352},
doi = {10.1002/mrm.26352},
pages = {142--155},
keywords = {dynamic magnetic resonance imaging, CMR, perfusion imaging, total generalized variation, infimal convolution, variational models},
year = {2017},
}
```


## Acknowledgement
This work is funded and supported by the [Austrian Science Fund (FWF)](http://fwf.ac.at) in the context of project 'SFB F32-N18' [Mathematical Optimization and Applications in Biomedical Sciences](http://imsc.uni-graz.at/mobis/)
and the province of Styria with the project 'HTI:Tech for Med (ABT08-22-T-7/2013-13)'. We also gratefully acknowledge the support of [NVIDIA Corporation](http://nvidia.com) with the donation 
of GPU computing hardware used for this research.

For questions and comments on the project please contact [Matthias Schloegl](mailto:matthias.schloegl@tugraz.at)
## Dependencies
* CUDA 4.0
* CMAKE 2.8
* GCC
* [AVIONIC](https://github.com/IMTtugraz/AVIONIC.git)
* [AGILE](https://github.com/IMTtugraz/AGILE.git)
* [gpuNUFFT](https://github.com/andyschwarzl/gpuNUFFT)
* [ISMRMRD](https://github.com/ismrmrd/ismrmrd)
* [DCMTK](http://dicom.offis.de/dcmtk.php.de)
* [Doxygen](http://www.stack.nl/~dimitri/doxygen/) (for code docs)
* [mapVBVD](https://github.com/cjohnevans/Gannet2.0/blob/master/mapVBVD.m) by Philipp Ehses (for reading Siemens raw data files, already included in repository)
* [MATLAB](https://www.mathworks.com/products/matlab.html) tested under the version R2016b

## Setup
0 Preparations
* build dcmtk from source (shared libs on) <br>
 `wget https://distfiles.macports.org/dcmtk/dcmtk-3.6.1_20160630.tar.gz`
* build hdf5 from source (shared libs on) <br>
 `wget https://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.10/src/hdf5-1.8.10.tar.bz2`
* Make sure that the CUDA environment is set up correctly 
* Ensure that the DCMDICTPATH environment variable is set correctly, namely with <br>
  `export DCMDICTPATH=$DCMDICTPATH:/usr/local/share/dcmtk/dicom.dic:/usr/local/share/dcmtk/private.dic`

1 Install AGILE lib 
```
git clone https://github.com/IMTtugraz/AGILE.git
cd AGILE
mkdir build
cd build
cmake ..
make -j 
sudo make install
``` 

2 Install gpuNUFFT 
```
git clone --branch deapo-scaling https://github.com/andyschwarzl/gpuNUFFT.git
cd gpuNUFFT/CUDA
mkdir build
cd build
cmake ..
make
``` 

3 Install ISMRMRD 
```
git clone https://github.com/ismrmrd/ismrmrd
cd ismrmrd/
mkdir build
cd build
cmake ../
make
sudo make install
``` 

4 Install AVIONIC recon lib
```
git clone https://github.com/IMTtugraz/AVIONIC.git
cd AVIONIC/CUDA
mkdir build
cd build
cmake .. -DGPUNUFFT_ROOT_DIR=/path/to/gpuNUFFT
make -j 
```

5 Add binary to PATH (bash)
```
in ~/.bashrc add:
export PATH=/path/to/AVIONIC/bin/:${PATH} 
```

6 Download BSReconFramework and download sample data
```
git clone https://github.com/IMTtugraz/BSReconFramework.git
cd ./BSReconFramework/data/
wget https://zenodo.org/record/1296051/files/gre_BlochSiegert_3D.dat
wget https://zenodo.org/record/1296051/files/gre_BlochSiegert_acc_12x4.dat
wget https://zenodo.org/record/1296051/files/smaps_walsh3d_slice.mat
```


## DEMO 1: Reconstruction of retrospectively subsampled data in the human brain using different subsampling patterns


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1296051.svg)](https://doi.org/10.5281/zenodo.1296051)

The file ``gre_BlochSiegert_3D.dat`` contains a fully-sampled Bloch Siegert 3D-dataset acquired with a GRE-sequence. For reconstruction run the file 

```
main_reconstructRetrospectivelySubsampled.m 
```

using MATLAB and select a subsampling pattern out of 

```
'full':     fully sampled
'block':    block pattern in k-space center
'eliptic':  eliptical pattern in k-space center
'vdrandom': variable density pattern
'gauss':    pattern with Gaussian density function
```

The resuls are written into a ``.mat`` file in the data folder as ``B1Map`` in µT and ``flipAngleMap`` as normalized nominal flip angle in %.


## DEMO 2: Reconstruction from a prospectively subsampled dataset in the human brain using a block pattern of size 12x4


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1296051.svg)](https://doi.org/10.5281/zenodo.1296051)

The file ``gre_BlochSiegert_acc_12x4.dat`` contains a prospectively subsampled Bloch Siegert 3D-dataset using a block pattern with size 12x4. For reconstruction run the file 

```
main_reconstructProspectivelySubsampled.m
```

using MATLAB. The resuls are written into a ``.mat`` file in the data folder as ``B1Map`` in µT and ``flipAngleMap`` as normalized nominal flip angle in %.


## Demo 3: Reconstruction of the fully sampled reference without any regularization


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1296051.svg)](https://doi.org/10.5281/zenodo.1296051)

For comparison the fully sampled reference can be gained out of ``gre_BlochSiegert_3D.dat`` by running the file

```
main_fullySampledRecon.m
```

using MATLAB. The resuls are written into a ``.mat`` file in the data folder as ``B1Map_full`` in µT and ``flipAngleMap_full`` as normalized nominal flip angle in %. 

For coil combination an implementation of the algorithm proposed by Walsh et al. was used. 
(Walsh DO, Gmitro AF, Marcellin MW. Adaptive reconstruction of phased array MR imagery. Magn Reson Med 2000;43(5):682–690) [DOI](https://doi.org/10.1002/(sici)1522-2594(200005)43:5<682::aid-mrm10>3.0.co;2-g)

