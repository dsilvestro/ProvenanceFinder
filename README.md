# ProvenanceFinder

### The methods implemented in ProvenanceFinder are described here:

Pérez-Consuegra, Nicolá., Parra, M., Jaramillo, C., Silvestro, D., Echeverri, Sebastiá., Montes, C., Jaramillo, José.Marí., Escobar, J., Provenance analysis of the Pliocene Ware Formation in the Guajira Peninsula, northern Colombia: Paleodrainage implications, Journal of South American Earth Sciences (2017), doi: 10.1016/j.jsames.2017.11.002.


### Probabilistic estimation of the most probable origin of detrital zircons among a set of potential sources.

ProvenanceFinder implements a probabilistic framework to infer the provenance of a set of zircons of unknown origin given a
set of potential source regions. A sample is defined here as the radiometric age of an individual zircon crystal. The term source refers to the mountain range from which a zircon could have been derived.
The main outputs of this analysis are as follows: 

(1) The likelihood of a set of samples (i.e., detrital zircon ages from a sedimentary unit, in this case the Ware Formation) given each source region. Likelihoods can be directly compared to assess which source best explains the distribution of ages measured in the analyzed set of samples.

(2) The relative probability for each sample to have originated in each of the possible sources. 

(3) The identification of samples that are not adequately explained by any of the samples in the range of sources. 

(4) The likelihood of a mixed model that allows each sample to have a different origin. This model does not assign each
sample to a source, but rather it integrates out the uncertainties around the provenance of the samples and
estimates the relative contribution of each source to the analyzed set of samples. 

All inferences about the age ranges explicitly incorporate the uncertainty around the dating of each sample (error for the individual age of a zircon crystal).

**ProvenanceFinder is licensed under a [AGPLv3 License](https://tldrlegal.com/license/gnu-affero-general-public-license-v3-(agpl-3.0)#summary).**

