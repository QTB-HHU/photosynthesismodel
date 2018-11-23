# Photosynthesis model
Computational model of photosynthesis resulted from a merge of two previously published mathematical models of relevant subprocesses. Original article available now as a pre-print: doi: https://doi.org/10.1101/476846 

## Model description
We have merged together two independent, previously developed kinetic models of photosynthesis, both based on ordinary differential equations. The first model describes the primary photosynthetic reactions through the photosynthetic electron transport chain (PETC), leading to the production of ATP and NADPH. It has been developed based on our previous work: the core model of the PETC by Ebenhöh *et al.*[1] and the model of high-energy dependent quenching in higher plants developed by Matuszyńska *et al.* [2]. The second model is the Poolman [3] implementation of the carbon fixation model by Pettersson and Ryde-Pettersson [4], reproduced in our Institute using the modelbase software [5]. 

## Code structure
The model has been implemented using modelbase package, a free expandable Python package for building and analysing dynamic mathematical models of biological systems [5]. The model is described in three files:
- [parameters.py](https://github.com/QTB-HHU/photosynthesismodel/blob/master/parameters.py), containing all parameters in the dot-access dict format [DotMap](https://github.com/drgrib/dotmap),
- [reactionrates.py](https://github.com/QTB-HHU/photosynthesismodel/blob/master/reactionrates.py), with functions defining the appropriate rate laws, and
- [model.py](https://github.com/QTB-HHU/photosynthesismodel/blob/master/model.py), containing information on model parameters, model variables, rate equations and stoichiometries.

The model object is instantiated by calling the
- [instantiate.py](https://github.com/QTB-HHU/photosynthesismodel/blob/master/instantiate.py).

## Repeat our analysis with the Jupyter Notebook
We are providing the Jupyter Notebook by Nima P. Saadat where all results included in the main manuscript can be easily rerun. We higly encourage you to start [with it](https://github.com/QTB-HHU/photosynthesismodel/blob/master/run.ipynb) while getting to know our model.

## References
[1] Ebenhöh O, Fucile G, Finazzi G, Rochaix JD, Goldschmidt-Clermont M. Short-term acclimation of the photosynthetic electron transfer chain to changing light: a mathematical model. *Philos Trans B* **2014**;369(1640):20130223. Open Access https://doi/10.1098/rstb.2013.0223 

[2] Matuszyńska A, Heidari S, Jahns P, Ebenhöh O. A mathematical model of non-photochemical quenching to study short-term light memory in plants. *Biochim Biophys Acta-Bioenerg* **2016**;1857(12):1–7. Open Access https://doi.org/10.1016/j.bbabio.2016.09.003

[3] Poolman MG, Fell DA, Thomas S. Modelling photosynthesis and its control. *J Exp Bot* **2000**;51(90001):319–328.

[4] Pettersson G,Ryde-Pettersson U. A mathematical model of the Calvin photosynthesis cycle. *Eur J Biochem* **1988**; 175(3):661–672.

[5] Ebenhöh O, Aalst vM, Saadat NP, Nies T, Matuszyńska A. Building mathematical models of biological systems with modelbase. *J Open Res Softw* **2018**;6(1):24. Open Access http://doi.org/10.5334/jors.236


