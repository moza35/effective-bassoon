# effective-bassoon
# This is a repository for a big data project I conducted in my second year of my BSc. I have prepared this for the LIDA data scientist role at the University of Leeds. Unfortunately, this is the only project that I previously worked on that I can share publicly since the others included human genomic data. 

The aim of this project was to investigate diversification of ancient crocodiles over time and its link to environmental change. For this project, I was provided with a macroevolutionary and environmental data set that explores the macroevolutionary drivers of a vertebrate clade called Pseudosuchia (crocodiles plus their extant & extinct relatives). This included a time-calibrated phylogenetic tree of Pseudosuchia (fossilCrocPhylogeny.tre).
diversification rate data for Pseudosuchia (fossilCrocDiversificationData.txt), and habitat data for Pseudosuchia (HabitatData.csv). It also included environmental data such as global temperature through time (temperatureTimeSeries.csv) and global sea level through time (seaLevelTimeSeries.csv).

Unfortunately, I can't upload this data to Git since the files are too large, but they're publicly available here: https://www-users.york.ac.uk/~kd856/WorkshopData/

For this project I plotted and analysed time series from environmental and diversification data using BAMM/fossilBAMM, conducted correlation analyses between time series, created phylogenetic trees using phytools, and created histograms and tested for statistical significance using ggplot2. 