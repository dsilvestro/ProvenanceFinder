# Import ProvenanceFinder functions
directory1 <- "/Users/myname/folder1/subfolder1/ProvenanceFinder.R"
source(directory1) 

#### INPUT FILE and MAIN SETTINGS
setwd("/Users/myname/folder1/subfolder1/") # set working directory

provenance_file  <- "provenance_example_file.txt"
name_target_sample  <- "unknown" # name of the sample with unknown origin. Same name 
                                 # as the one used in the provenance_example_file

#### ANALYSIS SETTINGS (default values)
n_points   <- 1000  # number of bins used to compute the probabilities 
                    # the higher the number, the more accurate 
                    # (and time consuming) the computation
plot_PDF   <- TRUE  # plot probability density funciton (PDF) of each source 
plot_CORR  <- TRUE  # plot correlations between PDFs of different sources
calcProb   <- TRUE  # calculate probability of target sample given each source (includes Mixed model)
runJK      <- 1000  # run Jackknife randomizations to assess robustness of the results (set to higher
	              # value for more reliable scores)

# run ProvenanceFinder analysis (fast run, change settings for more accurate results)
run_provenance_estimation(provenance_file,name_target_sample,n_points=1000,runJK=1000)

# The analysis produces 6 output files (save in the same directory as the input file):

"*_Summary.txt"               # main results: perobability of each source and mixed model, samples
                              # from unknown sources, and jackknife results
"*_prob_density_plots.pdf"    # probability density plot for each source (only if plot_PDF <- TRUE)
"*_logLikelihood_table.txt"   # log-likelihood of each sample conditional on each source 
"*_relProbability_table.txt"  # relative probability of each sample conditional on each source 
"*_source_correlation.pdf"    # plot correlations of probability densities between sources 
                              # (only if plot_CORR <- TRUE)
"*_prob_data.pdf"             # plot probability of the data (Equation 3)

