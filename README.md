# 9 Compartmental CO2 Model
This model and data adapts J. Gisolf et. al.'s work in their 2003 study ["Tidal volume, cardiac output and functional residual capacity determine end-tidal CO2 transient during standing up in humans"](https://pubmed.ncbi.nlm.nih.gov/14608002/). Data from the New Jersey Medical School Vestibular Autonomic Lab is used on this model. 

The goal of this project is to give a numeric approximation of the CO2 per lung compartment. In the process, the amount of end-tidal CO2 is estimated. The advantage of this method is that labs without the specific equipment to directly measure the CO2 per lung compartment can be estimated. 

The input data from the data include the estimated CO2. Many values in the model are taken directly from the paper, which may not apply precisely due to a different sample group.

Stroke Volume Per Breath had to be estimated using the [Liljestrand & Zander formula](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5317099/). Inaccuracies may be introduced this way. 

Feel free to contact me for more information about this project!