# functions for TR analysis

# PLOT THEMES FOR FLUDYNEMO FLU DATA ANLAYSIS:

SEGMENTS = c('H1N1_PB2','H1N1_PB1','H1N1_PA','H1N1_HA',
             'H1N1_NP','H1N1_NA','H1N1_MP','H1N1_NS')

 PlotTheme1 = theme_bw() +
            theme(legend.key = element_blank(),
                  strip.background = element_rect(colour="black", fill="white"),
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

PlotTheme2 = theme_bw() +
              theme(legend.key = element_blank(),
                   strip.background = element_rect(colour="black", fill="white"))

# completely blank backgrounds
PlotTheme3 = theme(legend.key = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"),
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                  strip.background = element_rect(colour="black", fill="white"))


###############################################################################
############################ FUNCTIONS #######################################
##################### Arrange Var (WITH rep, with stock AND ref info) ###########
prep_serial_metadata = function(metadata_df, remove_cals){

    # INPUT is the metadata df that was generated using the coverage notebook
    # there should be a column present that indicates whether or not the metadata was good
    # THIS WILL ONLY WORK FOR SERIAL SAMPLES!

    message("Number of samples pre filter: ", length(levels(factor(c(meta$name.y)))))

    # filter out bad samples. Filter out replicates where the library pools are the same
    # sequencing replicates were sequenced in different pool
    meta_filt = metadata_df %>% filter(quality == 'good' &
                           LIB_POOL.x != LIB_POOL.y) %>% droplevels()

    message("Number of samples after filtering for good samples: ", length(levels(factor(c(meta_filt$name.y)))))


    message("Remove Set 8 rep 2 comparisons")

    # redid the second replicate of set 8 samples in later seq runs due to contamination with rep2
    # selecting the pool 8. Replicate 2 had issues
    bad_s8 = meta_filt %>% filter(Set.x == 'S8' &
                                  Rep.x == 'R1' &
                                  Set.y == 'S8' &
                                  Rep.y == 'R2')

    # Pull out the names of the samples that were in set 8
    bad_s8_names = levels(factor(bad_s8$name.y)) # REMOVE THESE

    # filter the metadata for set 8 samples
    meta_filt = meta_filt %>% filter(!name.y %in% bad_s8_names) %>% droplevels()

    message("Number of samples post S8 filter: ", length(levels(factor(c(meta_filt$name.y)))))

    message("Determining which sample comparisons to use (p6 or p6 redo)")

    # the sequencer failed during the very end of the pool 6 run.
    # So pool 6 was repooled and resequenced
    # So here I am choosing which data to use (the better of the two runs)

    message("Determined that p6 should be used and NOT p6 redo")

    # filter first for pool p6 samples (only good samples here at this point)
    p6 = meta_filt %>% filter(LIB_POOL.x == 'P6' |
                            LIB_POOL.y == 'P6') %>% droplevels()

    # pull the ferret info
    p6_names = levels(factor(p6$ferret_info))

    # Filter for p6 redo samples (only good samples here at this point)
    p6r = meta_filt %>% filter(LIB_POOL.x == 'P6_REDO' | LIB_POOL.y == 'P6_REDO') %>% droplevels()

    # Pull names for the redos
    p6r_names = levels(factor(p6r$ferret_info))

    # Take the set difference of the two. Here, what comes out will indicate the redos were better
    p6r_better = setdiff(p6r_names, p6_names)

    # Take the set difference of the two. Here, what comes  out will indciate the first pool was better
    p6_better = setdiff(p6_names, p6r_names)

    # because p6 is all around better, we will just be using the p6 comparisons
    # filter out p6 redo samples:
    meta_filt = meta_filt %>% filter(LIB_POOL.x != 'P6_REDO' & LIB_POOL.y != "P6_REDO")

    meta_filt = meta_filt[!duplicated(meta_filt), ] %>% droplevels()

    message("Number of samples post p6 redo filter: ", length(levels(factor(c(meta_filt$name.y)))))

    message("Checking samples with more than one replicate: ")

    # group and pull out ferret information that is present in more than just two replicates (16 = 8 segments * 2 replicates)
    y = meta_filt %>% group_by_at(all_of(c("ferret_info"))) %>% tally() %>% filter(n>16)

    # filter for ^ those samples.
    y = meta_filt %>% filter(ferret_info %in% levels(factor(y$ferret_info)) &
                         Set.x == Set.y &
                         PlatePosition.x == PlatePosition.y &
                        Rep.x == 'R1' & Rep.y == 'R2')

    y = y[!duplicated(y), ] %>% droplevels()

    # going through MANUALLY and making a list of which to keep:
    keep_samples = c('S10_S10_C03_C03_P4_P5',
                     "S11_S11_C01_C01_P6_P7",'S11_S11_A01_A01_P6_P7',
                     'S10_S10_G07_G07_P4_P5', 'S10_S10_F03_F03_P4_P5',
                     'S11_S11_A05_A05_P6_P7', 'S10_S10_F09_F09_P4_P5',
                     'S10_S10_B08_B08_P4_P5', 'S11_S11_F01_F01_P6_P7',
                     'S17_S17_B04_B04_P6_P7', 'S11_S11_B03_B03_P6_P7',
                     'S16_S16_G12_G12_P6_P7', 'S10_S10_B06_B06_P4_P5',
                     'S10_S10_G05_G05_P4_P5',
                     'S8_S17_D04_E09_P4_P6',
                     'S10_S10_G05_G05_P4_P5')

    # removed 'S10_S10_B07_B07_P4_P5', ''S10_S10_A08_A08_P4_P5',' from kepit list after realizing it was a duplicate!!! 08.18.2021 duplicate

    y = y %>% filter(!name.y %in% keep_samples) # select the names that aren't in the keep vector

    remove_samps = levels(factor(y$name.y)) # these will be the samples names we want to remove

    meta_filt = meta_filt %>% filter(!name.y %in% remove_samps) %>% droplevels() # remove from the main metadata filter

    message("Number of samples post filter: ", length(levels(factor(c(meta_filt$name.y)))))

    z = meta_filt %>% group_by_at(all_of(c("ferret_info"))) %>% tally() %>% filter(n>16)

    z = meta_filt %>% filter(ferret_info %in% levels(factor(z$ferret_info)) & !name.y %in% keep_samples)

    remove_samps = levels(factor(z$name.y)) # remove the other samples replicates

    meta_filt = meta_filt %>% filter(!name.y %in% remove_samps) %>% droplevels()

    message("Number of samples post filter: ", length(levels(factor(c(meta_filt$name.y)))))

    # Are replicates both good? Are some bad?
    message("Check to make sure both replicates are good")

    only_one = meta_filt %>% group_by(name.y) %>% tally() %>% filter(n < 16)  # count the number of enteries for each sample. if both are good then there should be a total of 16 rows/name (8 segments, two replicates)

    message("Samples where both reps didn't pass thresholds: ", length(levels(factor(only_one$name.y))))

    also_remove = meta %>% filter(name.y %in% levels(factor(only_one$name.y)) & segment_count < 7) # filter for samples that we care about


    # NEED TO CHECK ON ?? tissue SAMPLE
    # FILTERING FOR JUST ONE STOCK SAMPLE SHOULD PROBABLY REARRANGE THIS!!!!

    meta_filt = meta_filt %>%
        filter(!name.y %in% c(levels(factor(also_remove$name.y))))%>%
        filter(tissue != '?') %>% filter(!name.y %in% remove_cals)

    message("Number of samples present after removing samples with bad reps: ", length(levels(factor(meta_filt$name.y))))

    return(meta_filt)
}

###############################################################################


ArrangeVarWRep = function(df){
    # remove replicate info because we do not have reps
    df = df %>% select('sample','segment','ntpos','nt',
                       'majmin','freq','aa','codon','nonsyn',
                       'binocheck','totalcount','aapos',
                       'refnt','stocknt','stockaa',
                       'freq1','freq2','totalcount1',
                       'totalcount2','nt1','nt2')

    # filter for major or minor info
    majvars = df %>% filter(majmin == 'major') %>% select(-nonsyn) %>% droplevels()

    minvars = df %>% filter(majmin == 'minor') %>% droplevels()

    # merge the information by columns that should be the same
    vdf = merge(majvars, minvars, by = c('sample','segment','ntpos',
                       'totalcount','aapos',
                       'refnt','stocknt','stockaa','binocheck',
                     'totalcount1','totalcount2') )

    # remove majmin columns
    vdf = vdf %>% select(-majmin.x, -majmin.y)

    # rename to be specific for major and minor information
    colnames(vdf) = c('sample','segment','ntpos','totalcount',
                  'aapos','refnt','stocknt','stockaa','binocheck',
                  'totalcount1', 'totalcount2',
                  'major','majorfreq','majoraa','majorcodon',
                  'majorfreq1','majorfreq2','major1','major2',
                  'minor','minorfreq',
                  'minoraa','minorcodon','nonsyn', 'minorfreq1',
                'minorfreq2','minor1','minor2')

    vdf$aapos = floor(vdf$aapos) # remove the 0.33, 0.66 from the aapos that was added when converting timo to python3

    vdf$aapos = as.numeric(as.character(vdf$aapos))

    return(vdf)
}


###############################################################################
##################### Arrange Var (no rep, with stock AND ref info) ###########

ArrangeVarNoRep = function(df){
    # remove replicate info because we do not have reps
    df = df %>% select('sample','segment','ntpos','nt',
                       'majmin','freq','aa','codon','nonsyn',
                       'binocheck','totalcount','aapos',
                       'refnt','stocknt','stockaa')

    # filter for major or minor info
    majvars = df %>% filter(majmin == 'major') %>% select(-nonsyn) %>% droplevels()

    minvars = df %>% filter(majmin == 'minor') %>% droplevels()

    # merge the information by columns that should be the same
    vdf = merge(majvars, minvars, by = c('sample','segment','ntpos',
                       'totalcount','aapos',
                       'refnt','stocknt','stockaa','binocheck') )

    # remove majmin columns
    vdf = vdf %>% select(-majmin.x, -majmin.y)

    # rename to be specific for major and minor information
    colnames(vdf) = c('sample','segment','ntpos','totalcount',
                  'aapos','refnt','stocknt','stockaa','binocheck',
                  'major','majorfreq','majoraa','majorcodon','minor','minorfreq',
                  'minoraa','minorcodon','nonsyn')

    return(vdf)

}


################################################################################
##################### Comparing cutoffs ######################################

CheckCutoffs = function(vdf){

    freqs = c(0.01, 0.02, 0.03, 0.05, 0.10)

    cuts = c(100, 200, 300, 500)

    fdf = data.frame()

    for (f in freqs){

        for (c in cuts){

            fdf_filt = FiltVar(vdf, minorfreq_cutoff = f,
               minor_coverage_cutoff = c,
               major_coverage_cutoff = c) # output is major AND minor. We will want to separate them

            message("Adding metadata that we care about")

            fdf_filt = AddMeta(fdf_filt, meta_filt %>%
                               select(name.y, ferret_id, day, tissue, cohort,
                                      exp_type, aged, sex, virus, study_type, model,
                                      study_id) %>% droplevels(), c('sample'), c('name.y'))

            fdf_filt = AddMeta(fdf_filt, sizes, c('segment'),c('segment'))

            message("Rearranging DF for figures")

            fdf_filt$segment = factor(fdf_filt$segment, levels = SEGMENTS)

            fdf_filt$tissue = factor(fdf_filt$tissue, levels = TISSUES)

            fdf_filt$day = factor(fdf_filt$day, levels = DAYLIST)

            fdf_filt$freq_cut = f

            fdf_filt$cov_cut = c

            fdf = rbind(fdf, fdf_filt)



            message("
                ")


        }

    }


    message("Separating into minor and maj dataframes and removing duplicates")

    min_df = fdf %>% filter(minorfreq >= min_freq & binocheck == 'True')

    min_df = min_df[!duplicated(min_df), ] %>% droplevels()


    message("Filtering min and maj dfs for NW, UL, LL, and CAL09 samples")

    keep_tissues = c('CAL09','NW','UL','LL')

    min_var = min_df %>% filter(minor %in% ntlist & major %in% ntlist) %>%
                         filter(tissue %in% keep_tissues &
                           aged == 'adult' &
                            day %in% DAYLIST | tissue == 'CAL09') %>% droplevels()


    return(min_var)

}


######### FUNCTIONS FOR METADATA AND COVERAGE ANALYSIS: ##################

##################### COVERAGE PLOTS #########################################
Covplot = function(df, savedir, coverage_cut, percentcov, prefix){

    xsamp = length(levels(factor(df$name)))

    #x = ggplot(data = df, aes(x=ntpos, y=totalcount, colour=(totalcount>=coverage_cut))) +
    x = ggplot(data = df, aes(x=ntpos, y=totalcount)) +

        geom_hline(yintercept = coverage_cut, linetype = 2) +

        geom_line(aes(group=1)) +

        #scale_color_manual(values=c('red','black')) +

        theme_bw() +

        theme(legend.key = element_blank(),
              strip.background = element_rect(colour="black", fill="white"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +

        #facet_grid(name~segment, scales ='free') +
        facet_grid(name~segment) +

        xlab("Nucleotide Position") +

        ylab("Raw Read Depth")

    #print(x)
    ggsave(x,
           filename = glue("{savedir}/{prefix}.CoveragePlot.{coverage_cut}.{percentcov}.pdf"),
           width = 12, height = xsamp*1.2, limitsize=FALSE)
}


# generate coverage files using the filtered mixed
Logcovplot = function(df, savedir, cov, percentcov, prefix){

    xsamp = length(levels(factor(df$name)))

    #x = ggplot(data = df, aes(x=ntpos, y=log10(totalcount), colour=(totalcount>=log10(cov)))) +
    x = ggplot(data = df, aes(x=ntpos, y=log10(totalcount))) +

        #geom_hline(yintercept = log10(cov), linetype = 2) +

        geom_hline(yintercept = log10(cov), linetype = 2) +

        geom_line(aes(group=1)) +

        #scale_color_manual(values=c('red','black')) +

        theme_bw() +

        theme(legend.key = element_blank(),
              strip.background = element_rect(colour="black", fill="white"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +

        #facet_grid(name~segment, scales ='free') +
        facet_grid(name~segment) +

        xlab("Nucleotide Position") +

        ylab("Raw Read Depth")

    #print(x)
    ggsave(x,
           filename = glue("{savedir}/{prefix}.CoveragePlot.log10.{cov}.{percentcov}.pdf"),
           width = 12, height = xsamp*1.2, limitsize=FALSE)
}

################################################################################

############################### CHECKING IF SAMPLES PASS THRESHOLDS##############
CoverageAcross = function(mydata, coverage_cut, percentcov, sizes, wkdir){
    # input: coverage file from pipeline, coverage cutoff, segment size information, expected library number
    # output: cov df that has the average coverage, percent cov above cutoff, and number of segs that pass cutoff for each sample
    # NEED A METADATA FP;DER AMD A COVERAGE_FIGURES FOLDER!
    # changed 06.23.2021 removed expected number of SAMPLES


    if (!dir.exists(glue("{wkdir}/coverage_figures"))) {
      dir.create(glue("{wkdir}/coverage_figures"))
      }

    if (!dir.exists(glue("{wkdir}/metadata"))) {
      dir.create(glue("{wkdir}/metadata"))
      }

    message(glue("Coverage cutoff is: {coverage_cut}x
                  Percentage covered cutoff is: {percentcov}%"))
    # calculate average coverage for segments across samples
    # this will help for samples that have low or 0 coverage
    Avg_it = aggregate(x = mydata$totalcount,
          by= list(name = mydata$name, segment= mydata$segment),
          FUN=mean)

    colnames(Avg_it) = c('name','segment','avg_coverage')

    # filter out positions that have coverage BELOW coverage cutoff
    # group by the id and segment
    # count the number of positions that are above the coverage cutoff
    Above_cut = mydata %>%
        filter(totalcount >= coverage_cut) %>%
        group_by(name, segment) %>%
        tally()

    mergeit = merge(Avg_it, Above_cut, by = c('name','segment'), all = TRUE)

    mergeit = mergeit[!duplicated(mergeit), ] %>% droplevels()

    mergeit$n[is.na(mergeit$n)] = 0  # if none above cut make it 0

    # merge with sizes so we can calculate how much of the genome is covered
    mergeit = merge(mergeit, sizes, by='segment',  all.x=T)

    # rename the columns
    colnames(mergeit)= c('segment', 'name', 'avg_coverage','above_cutoff','SegmentSize','STRAIN')

    # calculate the percentage of each segment that passes our coverage cutoff
    mergeit$percent_above = (mergeit$above_cutoff/mergeit$SegmentSize)*100

    # filter samples out those that have > 90% coverage in at least 1 segment
    per_df = mergeit %>% filter(percent_above >= percentcov) %>% droplevels()

    per_df = per_df[!duplicated(per_df), ] %>% droplevels() # remove any duplicates that may be there

    nopass = mergeit %>% filter(percent_above < percentcov) %>% droplevels()

    nopass = nopass[!duplicated(nopass), ] %>% droplevels()

    # count the number of segs for each sample that pass thresholds
    seg_count = per_df %>%
        group_by(name,STRAIN) %>%
    tally()

    colnames(seg_count) = c('name','STRAIN','segment_count')

    # Add the segment count information to other data, if there is an NA it means no segment passed the cutoffs
    mergeit = merge(mergeit, seg_count, by=c('name','STRAIN'), all = TRUE)

    mergeit$segment_count[is.na(mergeit$segment_count)] = 0

    # 'good' samples are those that have at least 200x coverage across 90% of all 8 segments
    good_samps = mergeit %>% filter(segment_count  == 8) %>% droplevels()

    good_samps = good_samps[!duplicated(good_samps), ] %>% droplevels()

    # good samps coverage plots

    Covplot(mydata %>% filter(name %in% levels(factor(good_samps$name))), glue("{wkdir}/coverage_figures/"), coverage_cut, percentcov, "GoodSamps")

    Logcovplot(mydata %>% filter(name %in% levels(factor(good_samps$name))), glue("{wkdir}/coverage_figures/"), coverage_cut, percentcov, "GoodSamps")

    message("Number of libraries that pass: ", nrow(good_samps)/8)

    bad = mergeit %>% filter( segment_count < 8 ) %>% droplevels()

    bad = bad[!duplicated(bad), ] %>% droplevels() # drop any duplicates that may come with rearranging dfs

    # bad samples plots here
    Covplot(mydata %>% filter(name %in% levels(factor(bad$name))), glue("{wkdir}/coverage_figures/"), coverage_cut, percentcov, "BadSamps")

    Logcovplot(mydata %>% filter(name %in% levels(factor(bad$name))), glue("{wkdir}/coverage_figures/"), coverage_cut, percentcov, "BadSamps")

    message("Number of libraries that DO NOT pass: ", nrow(bad)/8)

    message("Samples that did not pass coverage thresholds: ", list(levels(factor(bad$name))))

    good_samps$quality = 'good'

    bad$quality = 'bad'

    full_cov_df = rbind(good_samps, bad)

    full_cov_df = full_cov_df[!duplicated(full_cov_df), ] %>% droplevels()

    # check to make sure sample number is correct


    message("

    ")


    write.csv(full_cov_df,
         file=glue('{wkdir}/metadata/Coverage.meta.{coverage_cut}.{percentcov}.csv'),
         row.names=FALSE) ## the file name will have the path to save!

    return(full_cov_df)

}

CoverageAcrossControls = function(mydata, coverage_cut, percentcov, sizes, wkdir){
  # input: coverage file from pipeline, coverage cutoff, segment size information, expected library number
  # output: cov df that has the average coverage, percent cov above cutoff, and number of segs that pass cutoff for each sample
  # NEED A METADATA FP;DER AMD A COVERAGE_FIGURES FOLDER!
  # changed 06.23.2021 removed expected number of SAMPLES
  
  
  if (!dir.exists(glue("{wkdir}/coverage_figures"))) {
    dir.create(glue("{wkdir}/coverage_figures"))
  }
  
  if (!dir.exists(glue("{wkdir}/metadata"))) {
    dir.create(glue("{wkdir}/metadata"))
  }
  
  message(glue("Coverage cutoff is: {coverage_cut}x
                  Percentage covered cutoff is: {percentcov}%"))
  # calculate average coverage for segments across samples
  # this will help for samples that have low or 0 coverage
  Avg_it = aggregate(x = mydata$totalcount,
                     by= list(name = mydata$name, segment= mydata$segment),
                     FUN=mean)
  
  colnames(Avg_it) = c('name','segment','avg_coverage')
  
  # filter out positions that have coverage BELOW coverage cutoff
  # group by the id and segment
  # count the number of positions that are above the coverage cutoff
  Above_cut = mydata %>%
    filter(totalcount >= coverage_cut) %>%
    group_by(name, segment) %>%
    tally()
  
  mergeit = merge(Avg_it, Above_cut, by = c('name','segment'), all = TRUE)
  
  mergeit = mergeit[!duplicated(mergeit), ] %>% droplevels()
  
  mergeit$n[is.na(mergeit$n)] = 0  # if none above cut make it 0
  
  # merge with sizes so we can calculate how much of the genome is covered
  mergeit = merge(mergeit, sizes, by='segment',  all.x=T)
  
  # rename the columns
  colnames(mergeit)= c('segment', 'name', 'avg_coverage','above_cutoff','SegmentSize','STRAIN')
  
  # calculate the percentage of each segment that passes our coverage cutoff
  mergeit$percent_above = (mergeit$above_cutoff/mergeit$SegmentSize)*100
  
  # filter samples out those that have > 90% coverage in at least 1 segment
  per_df = mergeit %>% filter(percent_above >= percentcov) %>% droplevels()
  
  per_df = per_df[!duplicated(per_df), ] %>% droplevels() # remove any duplicates that may be there
  
  nopass = mergeit %>% filter(percent_above < percentcov) %>% droplevels()
  
  nopass = nopass[!duplicated(nopass), ] %>% droplevels()
  
  # count the number of segs for each sample that pass thresholds
  seg_count = per_df %>%
    group_by(name,STRAIN) %>%
    tally()
  
  colnames(seg_count) = c('name','STRAIN','segment_count')
  
  # Add the segment count information to other data, if there is an NA it means no segment passed the cutoffs
  mergeit = merge(mergeit, seg_count, by=c('name','STRAIN'), all = TRUE)
  
  mergeit$segment_count[is.na(mergeit$segment_count)] = 0
  
  # 'good' samples are those that have at least 200x coverage across 90% of all 8 segments
  good_samps = mergeit %>% filter(segment_count  == 3) %>% droplevels()
  
  good_samps = good_samps[!duplicated(good_samps), ] %>% droplevels()
  
  # good samps coverage plots
  
  Covplot(mydata %>% filter(name %in% levels(factor(good_samps$name))), glue("{wkdir}/coverage_figures/"), coverage_cut, percentcov, "GoodSamps")
  
  Logcovplot(mydata %>% filter(name %in% levels(factor(good_samps$name))), glue("{wkdir}/coverage_figures/"), coverage_cut, percentcov, "GoodSamps")
  
  message("Number of libraries that pass: ", nrow(good_samps)/3)
  
  bad = mergeit %>% filter( segment_count < 3 ) %>% droplevels()
  
  bad = bad[!duplicated(bad), ] %>% droplevels() # drop any duplicates that may come with rearranging dfs
  
  # bad samples plots here
  if(nrow(bad) > 0){
    Covplot(mydata %>% filter(name %in% levels(factor(bad$name))), glue("{wkdir}/coverage_figures/"), coverage_cut, percentcov, "BadSamps")
    
    Logcovplot(mydata %>% filter(name %in% levels(factor(bad$name))), glue("{wkdir}/coverage_figures/"), coverage_cut, percentcov, "BadSamps")
  
  } else{
    (message("Zero libaries DO NOT pass"))
  }
  
  message("Number of libraries that DO NOT pass: ", nrow(bad)/3)
  
  message("Samples that did not pass coverage thresholds: ", list(levels(factor(bad$name))))
  
  good_samps$quality = 'good'
  
  if(nrow(bad) > 0){
    bad$quality = 'bad'
  }
  
  full_cov_df = rbind(good_samps, bad)
  
  full_cov_df = full_cov_df[!duplicated(full_cov_df), ] %>% droplevels()
  
  # check to make sure sample number is correct
  
  
  message("

    ")
  
  
  write.csv(full_cov_df,
            file=glue('{wkdir}/metadata/Coverage.meta.{coverage_cut}.{percentcov}.csv'),
            row.names=FALSE) ## the file name will have the path to save!
  
  return(full_cov_df)
  
}


###############################################################################

FiltVar = function(df, minorfreq_cutoff = 0.03,
                   minor_coverage_cutoff = 200,
                  major_coverage_cutoff = 200){

    message("Frequency filter for minor variant: ", minorfreq_cutoff)
    message("Coverage cutoff for minor variant: ", minor_coverage_cutoff)
    message("Coverage cutoff for major variant: ", major_coverage_cutoff)

    # INPUT: REARRANGED DATAFRAME USING ARRANGEVARWREP or ARRANGEVAR FUNCTIONS
    # ONLY USEABLE WITH STOCK NT DATA!! NOT REF SPECIFIC

    # OUTPUT: FILTERED VARIANT DATAFRAME

    message("Input df length: ", nrow(df))

    message("Number of input samples: ", length(levels(factor(df$sample))))

    v = df %>% filter( minorfreq >= minorfreq_cutoff &
                              totalcount >= minor_coverage_cutoff &
                              binocheck == "True" |
                              stocknt != major &
                              totalcount >= major_coverage_cutoff)

    v = v[!duplicated(v), ] %>% droplevels()

    message("Output after filtering: ", nrow(v))

    message("Number of samples post filtering: ", length(levels(factor(v$sample))))

    return(v)

}

###############################################################################
################################# COUNT DF ####################################
TallyIt = function(df, groupit, new_colname){
    #INPUT: dataframe and vector of variables want to group by to count
    #OUTPUT: count dataframe using the group variables

    df = df[!duplicated(df), ] %>% droplevels

    count_df = df %>% group_by_at(vars(all_of(groupit))) %>% tally()

    colnames(count_df)[colnames(count_df) == 'n'] = new_colname

    return(count_df)
}

################################################################################

############### Adding metadata ###############################################

AddMeta = function(vcf_df, metadf, by_vcf, by_meta){

    temp = merge(vcf_df, metadf, by.x=all_of(by_vcf), by.y=all_of(by_meta), all.x = TRUE)

    temp = temp[!duplicated(temp), ] %>% droplevels()

    return(temp)

}

################################################################################
#05.10.2020: added a keep samps function for stock samples!
DistanceMatrices_keep = function(distfiles, distmethod, wkdir, goodsamps, keepsamps){

   for (filename in distfiles) {

       print(filename)

       mydata = read.csv(file=filename,header=T,sep=",",na.strings = c(''), row.names = NULL) # read in data

       mydata_filt = mydata %>% filter(sample %in% goodsamps) %>%
               select(sample,ntpos,A,C,G,T) %>%
               droplevels() # remove bad samples with low coverage, and select nt data

       mydata_filt = mydata %>% filter(sample %in% keepsamps) %>%
                       select(sample,ntpos,A,C,G,T) %>%
                       droplevels() # remove bad samples with low coverage, and select nt data

       set2 = unique(factor(mydata_filt$sample)) # used to check # of samples

       print(length(set2))

       dist_orig = matrix(data = 0, nrow = length(set2) , ncol = length(set2)) # generate empty matrix to add dist numbers to

       samples = unique(sort(mydata_filt$sample)) # sort so samples always in order

       rownames(dist_orig) = samples

       colnames(dist_orig) = samples

       variant_positions = c(unique(mydata_filt$ntpos)) # generate ntpos list

       message("Size of segment matrix to add to: ", dim(dist_orig))

       for (nt in variant_positions) {

           # iterate through ntpos
           d1 = (data=subset(mydata_filt, ntpos == nt, c("sample","A","C","G","T")))

           d1 = d1[!duplicated(d1),] %>% droplevels() # remove any dups

           set1 = unique(factor(d1$sample)) #check to make sure everything is the same size throughout

           setdf = setdiff(set2,set1) #if it isn't the same print what is different between the files

           #to make first column with samples names the row names
           rownames(d1) = d1$sample # make rownames the sample names

           d1[,1] = NULL # remove the row names columns (which is the first column selected)

           # [ row, columns]
           d1 = d1[order(row.names(d1)),] # order the rows

           #calculate distance
           #can change method- several distance calculations, try 'manhattan'/L1
           D_euc = dist(d1, method=distmethod)

           D_euc = as.matrix(D_euc) #change dist type to matrix

           # check to make sure dist_orig and dist_euc are in same order
           # change order of dist_orig to order of dist_euc
           dist_orig = dist_orig[rownames(D_euc), ]

           dist_orig = dist_orig[,colnames(D_euc)]

           #print true or false if they are in the same order
           #message("row order: ", all(rownames(dist_orig) == rownames(D_euc)))
           #message("row order of do and colnames of deu: ", all(rownames(dist_orig) == colnames(D_euc)))
           #message("row order of deu and colnames of do: ", all(colnames(dist_orig) == rownames(D_euc)))


           dist_orig <- dist_orig + D_euc #add to the overall matrix
         }

       dist.df = as.data.frame(dist_orig)

       dist.df$name = rownames(dist.df)

       seg_dist = dist.df %>% filter(is.na(dist.df) != 'TRUE')

       write.csv(dist.df,
            file=glue('{wkdir}{filename}.matrix.csv'),
            row.names=FALSE) ## the file name will have the path to save!

   }
}


DistanceMatrices = function(distfiles, distmethod, wkdir, goodsamps){

  # updated 06.28.2021 to filter FOR good
  # updated 06.28.2021: changed matrix addition distmethod

  message("Distance matrices")

  bad_samps = c()

   for (filename in distfiles) {

       message("Input file: ", filename)

       mydata = read.csv(file=filename,header=T,sep=",",na.strings = c(''), row.names = NULL) # read in data

       mydata_filt = mydata %>% filter(sample %in% goodsamps) %>%
               select(sample,ntpos,A,C,G,T) %>%
               droplevels() # remove bad samples with low coverage, and select nt data

       set2 = unique(factor(mydata_filt$sample)) # used to check # of samples

       message("Number of input samples: ", length(set2))

       dist_orig = matrix()

       variant_positions = c(unique(mydata_filt$ntpos)) # generate ntpos list

       message("Length of segment: ", length(variant_positions))

       for (nt in variant_positions) {

           d1 = (data=subset(mydata_filt, ntpos == nt, c("sample","A","C","G","T")))  # filter by ntpos

           d1 = d1[!duplicated(d1),] %>% droplevels() # remove any dups

           set1 = unique(factor(d1$sample)) #check to make sure everything is the same size throughout

           setdf = setdiff(set2,set1) #if it isn't the same print what is different between the files

           #to make first column with samples names the row names
           rownames(d1) = d1$sample # make rownames the sample names

           d1 = d1 %>% select(-sample)

           #can change method- several distance calculations, try 'manhattan'/L1
           D_euc = dist(d1, method=distmethod)

           D_euc = as.matrix(D_euc) #change dist type to matrix

           if (all(is.na(dist_orig))){

             dist_orig = matrix(data = 0, nrow = length(set2) , ncol = length(set2))

             rownames(dist_orig) = rownames(D_euc)

             colnames(dist_orig) = colnames(D_euc)

             dist_orig[rownames(D_euc), colnames(D_euc)] =  dist_orig[rownames(D_euc), colnames(D_euc)] + D_euc

              }else(dist_orig[rownames(D_euc), colnames(D_euc)] = dist_orig[rownames(D_euc), colnames(D_euc)] + D_euc)

       }
    dist.df = as.data.frame(dist_orig)

    dist.df$name = rownames(dist.df)

    remove_samps = levels(factor(rownames(dist.df %>% filter(is.na(dist.df) == 'TRUE'))))

    bad_samps = c(bad_samps, remove_samps)

    #dist.df = dist.df %>% filter(!name %in% remove_samps) %>% select(-(all_of(remove_samps)))

    message("The following samples were REMOVED from the matrix. Either too low of coverage OR there was a mismatch in the replicate files.You need to go back and check to determine.:
    ", list(remove_samps))

    message("Saving as matrix file")

    write.csv(dist.df,
            file=glue('{wkdir}/{filename}.matrix.csv'),
            row.names=FALSE) ## the file name will have the path to save!

    }
  return(levels(factor(bad_samps)))
}


#############################################################################
######################## FLUDYNEMO STOCKS MDS PLOTS!! ###############################
# SPECIFIC TO THIS PROJECT ONLY!
SegMDSPlots = function(df,filename){
   # SPECIFIC FOR TR FLUDYNEMO TISSUE DATA!


   p1 = ggplot(df, aes(x= mds1, y = mds2, color=day, shape = tissue)) +
           geom_point() +
           ggtitle(filename) +
           theme_bw() +
        PlotTheme1

   print(p1)

   return(p1)


}
##############################################################################
################# SUM DISTANCE MATRICES #####################################

#  06.28.2021 updated summing

SumMats = function(distmats, wkdir, metadata, metaby, badsamps){

   dist_full = matrix() # generate empty matrix to add dist numbers to

   for (matfile in distmats){

       print(matfile)

       mydata = read.csv(file=matfile,header=T,sep=",",na.strings = c('')) # read in data

       mydata = mydata %>% filter(!name %in% badsamps) %>% select(-all_of(badsamps)) %>% droplevels()

       rownames(mydata) = mydata$name  ## add rownames from name column

       sampnumber = length(levels(factor(mydata$name)))

       print(length(levels(factor(mydata$name))))

       mydata = mydata %>% select(-name)  ## remove name column

       mydata = apply(as.matrix(mydata[ , ]), 2, as.numeric)  # change to matrix make numeric values even wif na is present

       rownames(mydata) = colnames(mydata)

       print(dim(mydata))

       #mydata = data.matrix(mydata) ## turn into matrix

       if (all(is.na(dist_full))){

         dist_full = matrix(data = 0, nrow = sampnumber , ncol = sampnumber)

         rownames(dist_full) = rownames(mydata)

         colnames(dist_full) = colnames(mydata)

         dist_full[rownames(mydata), colnames(mydata)] =  dist_full[rownames(mydata), colnames(mydata)] + mydata


       }else(dist_full[rownames(mydata), colnames(mydata)] = dist_full[rownames(mydata), colnames(mydata)] + mydata)



       #dist_full = dist_full + mydata

       fit = cmdscale(mydata, eig = TRUE, k=2)

       mds1 = fit$points[,1]

       mds2 = fit$points[,2]

       sample = row.names(mydata)

       df = cbind(sample, mds1)

       df = cbind(df, mds2)

       df = as.data.frame(df)

       df$mds1 = as.numeric(as.character(df$mds1))

       df$mds2 = as.numeric(as.character(df$mds2))

       df = merge(df, metadata, by.x=c('sample'),
                by.y=c(metaby)) %>%
               droplevels()

       df = df[!duplicated(df), ] %>% droplevels()

       #print(colnames(df))

       SegMDSPlots(df, matfile)

}

 return(dist_full)

 fit = cmdscale(mydata, eig = TRUE, k=2)

 mds1 = fit$points[,1]

 mds2 = fit$points[,2]

 sample = row.names(mydata)

 df = cbind(sample, mds1)

 df = cbind(df, mds2)

 df = as.data.frame(df)

 df$mds1 = as.numeric(as.character(df$mds1))

 df$mds2 = as.numeric(as.character(df$mds2))

 df = merge(df, metadata, by.x=c('sample'),
        by.y=c(metaby)) %>%
       droplevels()


 df = df[!duplicated(df), ] %>% droplevels()

 dist_full = as.data.frame(dist_full)

 dist_full$name = rownames(dist_full)

 write.csv(dist_full,
         glue("{wkdir}distance_mats/H1N1.fullgenome.distance.matfile.csv"), row.names = FALSE)

}
###############################################################################
############################# MDS OF FULL DIST FILE ##########################
MDSfull = function(mydata, metadata, metaby){

   fit = cmdscale(mydata, eig = TRUE, k=2)

   mds1 = fit$points[,1]

   mds2 = fit$points[,2]

   sample = row.names(mydata)

   df = cbind(sample, mds1)

   df = cbind(df, mds2)

   df = as.data.frame(df)

   df$mds1 = as.numeric(as.character(df$mds1))

   df$mds2 = as.numeric(as.character(df$mds2))

   df = merge(df, metadata, by.x=c('sample'),
        by.y=c(metaby)) %>%
       droplevels()

   df = df[!duplicated(df), ] %>% droplevels()

   #print(head(df))

   return(df)

}


#############################################################################














































































################################FERRET COUNTS #################################

CheckNumbers = function(df, prefix, savedir){
  # INPUT: variant dataframe
  # OUTPUT: figures with counts of ferrets per day per facet vars

    df_filt = df %>% select(FERRET_ID, DAY, PRECHALLENGE, PRE_EXPOSURE, GENERATION) %>% droplevels()

    df_filt = df_filt[!duplicated(df_filt), ] %>% droplevels()

    df_filt$DAY = factor(df_filt$DAY, levels = daylist1)

    cplot = ggplot(df_filt, aes(x=DAY)) +
        geom_bar(color = 'black') +
        ylab("Number of ferrets") +
        xlab("DPI") +
        geom_hline(yintercept = 3, linetype = 2, color = 'red') +
        facet_grid(.~PRE_EXPOSURE + GENERATION, scales = 'free_x',space='free') +
        theme_bw() +
        theme(legend.key = element_blank(),
              strip.background = element_rect(colour="black", fill="white"))

    print(cplot)

    ggsave(cplot,
       filename = glue("{savedir}/{prefix}.ByPreExposure.pdf"),
       width = 9,
       height = 5, limitsize=FALSE)

    cplot2 = ggplot(df_filt, aes(x=DAY)) +
        geom_bar(color = 'black') +
        ylab("Number of ferrets") +
        xlab("DPI") +
        geom_hline(yintercept = 3, linetype = 2, color = 'red') +
        facet_grid(.~PRECHALLENGE + GENERATION, scales = 'free_x',space='free') +
        theme_bw() +
        theme(legend.key = element_blank(),
              strip.background = element_rect(colour="black", fill="white"))

    ggsave(cplot2,
       filename = glue("{savedir}/{prefix}.Prechallenge.pdf"),
       width = 6,
       height = 5, limitsize=FALSE)

    print(cplot2)
}

################################################################################

############################### SNV ANALYSIS ##################################


############################## TAKE MEAN, SD, SE, COUNT, MIN, MAX OF DF ########
MeanIt = function(df, groupitby, mean_col){
    library(plotrix)

    df = df[!duplicated(df), ] %>% droplevels()

    colidx = which(colnames(df)==mean_col)

    mean_df =  df %>%
        group_by_at(all_of(groupitby)) %>%
        mutate(count = n()) %>%
        group_by_at(all_of(c(groupitby, "count"))) %>%
        summarise_at(.vars = colnames(.)[colidx],
                   .funs = lst(min, max, mean = mean, median = median,
                              sd = sd, se = std.error))
    return(mean_df)
}

###############################################################################

###################### KRUSKAL WALLIS RANK SUM TEST FOR SIG ####################
KWrsTest = function(df, compcol, groupcol, pvalcut){

  colidx = which(colnames(df)==groupcol)

  colidx2 = which(colnames(df)==compcol)

  g = list(levels(factor(df[[colidx]])))

  kpval = kruskal.test(df[[colidx2]] ~ df[[colidx]], data = df)$p.value

  if (kpval <= pvalcut){
      message(glue("Significant difference between the groups (Kruskal-Wallis, pcutoff used: {pvalcut})

  pval: {kpval}

  groups: {g}"))
  } else{
      message(glue("*NO* significant difference between the groups (Kruskal-Wallis, pcutoff used: {pvalcut})

  pval: {kpval}

  groups: {g}"))
  }

}

############################################################################################################

##################### PAIRWISE WILCOX TEST ##################################
PairWCoxtestBon = function(df, compcol, groupcol){
  colidx = which(colnames(df)==groupcol)

  colidx2 = which(colnames(df)==compcol)

  g = list(levels(factor(df[[colidx]])))

  message(glue("formula for pairwise wilcox test is: {compcol} ~ {groupcol}"))

  message("")

  message("Vars used for test: ", list(g))

  pairwise.wilcox.test(df[[colidx2]], df[[colidx]],
                 p.adjust.method = "bonferroni", correct = FALSE)

}

################################################################################

################ PAIRWISE T-TEST #############################################
PairTtestBon = function(df, compcol, groupcol){
  colidx = which(colnames(df)==groupcol)

  colidx2 = which(colnames(df)==compcol)

  g = list(levels(factor(df[[colidx]])))

  message(glue("formula for pairwise t-test is: {compcol} ~ {groupcol}"))

  message("")

  message("Vars used for test: ", g)

  pairwise.t.test(df[[colidx2]], df[[colidx]], p.adj = "bonf")

}

###############################################################################
####################### SHANNON AT POSITION ###################################
ShannonPosMajmin = function(df){

  ## uses unformated csv variant files
  ## this will also use positions where there wasn't a minor variant - misleading

  df$majmin_shan = -(df$freq)*(log2(df$freq))

  df = df %>% group_by(sample, segment, ntpos) %>% mutate(shannon_ntpos = sum(majmin_shan))

  df = df %>% group_by(sample, segment) %>% mutate(segment_shan = sum(shannon_ntpos))

  df = df %>% group_by(sample) %>% mutate(shannon = sum(segment_shan))

  return(df)
}

ShannonPos = function(df){
  ## uses reformated variant info
  ## only uses positions with minorvariants
  df$shannon_ntpos = (-(df$majorfreq)*(log2(df$majorfreq))) + (-(df$minorfreq)*(log2(df$minorfreq)))

  df = df %>% group_by(sample, segment) %>% mutate(segment_shan = sum(shannon_ntpos))

  df = df %>% group_by(sample) %>% mutate(shannon = sum(segment_shan))

  return(df)
}
###############################################################################
########################### MELT DISTANCE ###################################

MeltDistance = function(dist_df, meta, metacol){

    meta = meta[!duplicated(meta), ] %>% droplevels()

    dist_melt = melt(dist_df)

    dist_melt = merge(dist_melt, meta, by.x = c('Var1'), by.y = metacol,
                     all.x = TRUE) %>% droplevels()

    dist_melt = merge(dist_melt, meta, by.x=c("Var2"), by.y = metacol,
                     all.x = TRUE)

    dist_melt = dist_melt[!duplicated(dist_melt), ] %>% droplevels()

    return(dist_melt)
}



###############################################################################
############################## DVG ANALYSIS #################################

DVG_FiltArrange = function(df, cov, count_cut=10){

  names(cov)[names(cov) == 'totalcount'] = 'covcount'

  message("DVG filter is set at: ", count_cut)

  message("Number of rows before filter: ", nrow(df))

  df = df %>% filter(DVG_freq >= count_cut) %>% droplevels()

  message(glue("Number of rows post filter of {count_cut}: "), nrow(df))

  df$bound = df$GroupBoundaries

  df = df %>% separate(bound, c("s1","s2","e1","e2"), sep = "([-_])")

  df = df %>% select(name, segment, segment_size, totalreads, DVG_group,
                   NewGap, NewStart, NewEnd, GroupBoundaries,
                   DeletionSize, DVG_freq, s1, s2, e1, e2) %>% droplevels()

  df_updated = df %>% group_by(name, segment) %>% mutate(SegmentDVGCount = sum(DVG_freq))

  message("Removing prefiltered df from memory")

  rm(df)
  # calculate percentage of diff dvgs to total dvgs
  df_updated$DVG_PercentTotalDVG = (df_updated$DVG_freq/df_updated$SegmentDVGCount)*100

  # caclulate percent abund out of total fragments aligned:
  df_updated$DVG_PercentAbundance = (df_updated$DVG_freq/df_updated$totalreads)*100

  df_updated$s1_1 = as.numeric(as.character(df_updated$s1)) - 1

  df_updated$e2_1 = as.numeric(as.character(df_updated$e2)) + 1

  message("Merging DVG df with coverage information")

  df_updated = merge(df_updated, cov, by.x=c('name','segment','s1_1'),
                 by.y = c("name","segment","ntpos"), all.x = TRUE)

  names(df_updated)[names(df_updated) == 'covcount'] = 'start_minus1'

  df_updated = merge(df_updated, cov, by.x=c('name','segment','e2_1'),
                 by.y = c("name","segment","ntpos"), all.x = TRUE)

  names(df_updated)[names(df_updated) == 'covcount'] = 'end_plus1'

  df_updated = merge(df_updated, cov, by.x=c('name','segment','s1'),
                 by.y = c("name","segment","ntpos"), all.x = TRUE)

  names(df_updated)[names(df_updated) == 'covcount'] = 'start_count'

  df_updated = merge(df_updated, cov, by.x=c('name','segment','e2'),
                 by.y = c("name","segment","ntpos"), all.x = TRUE)

  names(df_updated)[names(df_updated) == 'covcount'] = 'end_count'

  message("Making final calculations")

  df_updated$avg_count = (df_updated$start_minus1 + df_updated$end_plus1 + df_updated$start_count + df_updated$end_count)/4

  # using coverage to calculate frequency information
  df_updated$DVG_position_freq = (df_updated$DVG_freq/df_updated$avg_count)*100

  # dvg fpkm
  df_updated$DVG_FPKM = df_updated$DVG_freq/((df_updated$segment_size/1000)/(df_updated$totalreads/1000000))

  # calculate percent abundance out of total aligned to segment
  df_updated$TotDVGPercentAbund = (df_updated$SegmentDVGCount/(df_updated$totalreads/(df_updated$segment_size/1000)))*100

  df_updated$Segment_FPKM = df_updated$SegmentDVGCount/((df_updated$segment_size/1000)/(df_updated$totalreads/1000000))

  message("Removing coverage df from memory")

  rm(cov)

  return(df_updated)
}

##############################################################################
###############################################################################
############# ADD DVG MIDPOINT #############################################


MidPoint = function(df){

  message("Adding DVG midpoints")

  df$GapMidPoint = ((df$NewEnd - df$NewStart)/2) + df$NewStart

  df$DistanceFromSegCent =  df$GapMidPoint - (df$segment_size / 2)

  df$SegmentCenter = df$segment_size/2

  return(df)
}

##############################################################################
################### AVG DVG REPS ########################################
AvgDVGReps = function(df){

    head(df)
    df$DVG_freq = (df$DVG_freq.x + df$DVG_freq.y)/2

    df$DVG_PercentTotalDVG = (df$DVG_PercentTotalDVG.x + df$DVG_PercentTotalDVG.y)/2

    df$DVG_PercentAbundance = (df$DVG_PercentAbundance.x + df$DVG_PercentAbundance.y)/2

    df$avg_count = (df$avg_count.x + df$avg_count.y)/2

    df$DVG_position_freq = (df$DVG_position_freq.x + df$DVG_position_freq.y)/2

    df$DVG_FPKM = (df$DVG_FPKM.x + df$DVG_FPKM.y)/2

    df$TotDVGPercentAbund = (df$TotDVGPercentAbund.x + df$TotDVGPercentAbund.y)/2

    df$Segment_FPKM = (df$Segment_FPKM.x + df$Segment_FPKM.y)/2

    return(df)

    }
################################################################################


##############################################################################

flattenCorrMatrix <- function(cormat, pmat) {
    # taken from: http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
    ut <- upper.tri(cormat)
    data.frame(
        row = rownames(cormat)[row(cormat)[ut]],
        column = rownames(cormat)[col(cormat)[ut]],
        cor  =(cormat)[ut],
        p = pmat[ut]
    )
}

##############################################################################
