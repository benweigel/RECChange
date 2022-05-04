**RECChange**

This repository accompanies the manuscript **"Climate change reshuffles northern species within their niches"** by Antao & Weigel et al. and provides all necessary R scripts and data to reproduce included results.

**Antão, L.H.** * **& Weigel, B.** *, Strona, G., Hällfors, M., Kaarlejärvi, E., Dallas, T., Opedal, Ø., Heliölä, J., Henttonen, H., Huitu, O., Korpimäki, E., Kuussaari, M., Lehikoinen, A., Leinonen, R., Lindén, A., Merilä, P., Pietiäinen, H., Pöyry, J., Salemaa, M., Tonteri, T., Vuorio, K., Ovaskainen, O., Saastamoinen, M., Vanhatalo, J., Roslin, T., Laine, A-L. Climate change reshuffles northern species within their niches. *Nature Climate Change* (2022).


We used four decades of occurrence data on 1,478 species of birds, mammals, small rodents, butterflies, moths, forest understory vascular plants, and freshwater phytoplankton from Finland, to analyse shifts in species’ relative niche position over time – i.e., whether species occur at the lower end, at the optimum, or at the upper end of their climatic niche (each dataset was collected and treated independently). In addition, we analysed whether the relative importance of different climatic variables varied over time. As niche dimensions, we included annual values of mean temperature, total precipitation, duration of snow cover, and the North Atlantic Oscillation (NAO) index. We fitted joint species distribution models for each combination of taxonomic group × bioclimatic zone × decade in our data.

**Disclaimer**: The code in this repository represents one version of the code developed for the project and may yet undergo changes and revisions.

**Authors**: The code was developed through collaboration between Benjamin Weigel, Laura Antao, Giovanni Strona, Jarno Vanhatalo and Otso Ovaskainen.


**Contacts:**     
Benjamin Weigel - benjamin.weigel@helsinki.fi     
Laura Antao - laura.antao@helsinki.fi


**Data Files**

*Species data*    
**Bird**  ...   
**Butterfly** ...   
**Moth** monitoring data is available through the Finnish Biodiversity Information Facility, FinBIF (https://laji.fi/en/observation/list?sourceId=KE.1501)    
**Phytoplankton** data was obtained from the open access data service of Finnish Environment Institute (SYKE) (Open data portal http://www.syke.fi/en-US/Open_information)    
For the data sets **mammals, small rodents,** and **understory plants** belonging to the Natural Resources Institute Finland (Luke), any data requests should be sent to kirjaamo@luke.fi.    Raw data is available from data owners on reasonable request.

*Climatic data*  
We extracted values of daily mean temperature, daily precipitation sum, and daily snow depth from the Finnish
Meteorological Institute(https://etsin.fairdata.fi/datasets/fmi?keys=Finnish%20Meteorological%20Insitute&terms=organization_name_en.keyword&p=1&sort=best; first accessed in April 2019 and updated in May 2020). NAO values were derived from: NAO Index Data provided by the Climate Analysis Section, NCAR, Boulder, USA, Hurrell (2003). Updated regularly. Accessed 16.12.2019. https://climatedataguide.ucar.edu/climate-data/hurrell-north-atlantic-oscillation-nao-index-pc-based


Raw data on species monitoring can be accessed via the links in the manuscript, and/or upon reasonable request to the data owners.
The data shared here are the processed data used in the analysis, that is occurrence records after the data inclusion criteria described in the Methods.
For mammals, small rodents, and understory plant datasets, the accuracy of the coordinates for the sampling sites was coarsened to comply with data owners’ requirements. The files contain a “Sample ID” variable that identifies each sample, i.e. species occurrence for each site in each year.
The files also contain the climatic variables extracted for each site and year (annual mean temperature, annual summed precipitation, annual snow cover duration days, and NAO).  

The **data folder** contains folders for each taxonomic group. Each taxon folder includes .csv files with all data to fit each model, i.e. bioclimatic zone x decade combination. Csv files are labled after taxon_decade_zone.csv. Script *01 HMSC_RECC.R* serves as illustrative example script to fit the model for birds in the south boreal zone in decade 1, using *birds_decade_1_SB.csv* as data input. When wanting to iterate through several/all input files, one needs to adapt the code and make sure to rename the model file/output accordingly.


**R Analysis Files**    
The scripts exemplify how to fit the models for each taxon in each zone and decade (01), score species position within niche domains (02), extract the relevant parameter estimates from the JSDMs (03), and calculate the species composition (beta-diversity) metrics (04).

Please note that some of the code in this repository was written to run on a HPC cluster and that fitting the JSDMs can take multiple weeks for the largest subsets.

**Requirements**    
Program R version used 3.6.3 (2020-02-29)   
Other attached packages: *betapart_1.5.4, forcats_0.5.0, stringr_1.4.0, dplyr_1.0.7, purrr_0.3.4, readr_1.4.0, tidyr_1.1.4, tibble_3.1.6, ggplot2_3.3.5, tidyverse_1.3.0, Hmsc_3.0-12, coda_0.19-4*
