# InverseCardiacOutOfSample

## Table of contents
* [General info](#general-info)
* [Technologies](#technologies)
* [Code Examples](#code-examples)
* [Acknowledgement](#acknowledgement)

## General info
This project contains the prototype codes for computing the out-of-sample tuning of regularization parameter on inverse cardiac problem. These codes reproduce parts of the experiments reported in following paper:
* F. M. Melgarejo-Meseguer, E. Everss-Villalba, M. Gutierrez-Fernandez-Calvillo, S. Munoz-Romero, F.-J. Gimeno-Blanes, A. Garcia-Alberola, and J.-L. Rojo-Alvarez, “Generalization and Regularization for Inverse Cardiac Estimators,” IEEE Transactions on Biomedical Engineering, vol. 69, no. 10. Institute of Electrical and Electronics Engineers (IEEE), pp. 3029–3038, DOI:10.1109/TBME.2022.3159733
	
## Technologies
Project was created with:
* MATALAB: 2021.b


## Code Examples
The main file from the code is called "ToyExampelOutOfSample" and it allows us to reproduce two different experiments. The first one performs the out of sample regularization tuning on a simple 2D substrate with a Gaussian pulse acting as stimulation. After running this part of the code, a figure divided in 6 subpanels is generated, the panels show from left to right and up to down the measured extracellular potential, the estimed one, the difference among measured and estimed extracellular potential, the measured transmembrane potential, the estimed one, and the difference among measured and estimed transmembrane potential. This part of code goes up to line 121.

The second part of the code, which goes from line 122 to the end, performs the out of sample regularization tuning on an example from EDGAR database. After running this part of the code, two different figures are generated. The first one, which is divided in 4 different panels, shows on first row panels the Person's correlation coefficient and the RMSE computed among the measured and estimed cardiac potentials, and second row panels exhibit the same computations among the measured and estimed torso potentials. The second figure is divided in 6 different panels, they show in the first row the real and estimed signals for the worst correlation coefficient, the worst RMSE value, and the best correlation coefficient among measured and estimed torso potentials, and in the second row they show the same case but computed for cardiac potentials.

## Acknowledgement

* SCI Institute, “SCIRun: A. Scientific computing problem solving envi-
ronment, scientific computing and imaging institute (SCI),” 2016. [On-
line]. Available: http://www.scirun.org

* J. Tate et al., “Reducing error in ECG forward simulations with improved
source sampling,” Front. Physiol., vol. 9, 2018, Art. no. 1304.
