browser_plotter <-
    function(
        regions.gr,
        dataset.grlist,
        binsize = 100,
        bin_FUN = sum,
        pad_left = 1000,
        pad_right = 1000,
        scale_bar_size = 5000,
        color_plus = "#BB0021",
        color_minus = "#3B4992",
        color_unstranded = "black",
        tx_name = "tx_name",
        .labs = as.character(unlist(lapply(dataset.grlist, names))),
        .expand_ranges = setNames(rep(TRUE, length(names(dataset.grlist))), names(dataset.grlist))
    ) {
        # This function generates a 'browser shot' from a GRanges that includes a
        # single gene, and a list of lists where each sublist is a list of GRanges
        # of a consistant type of signal (stranded/unstranded).
        # For example, the structure could be a list of 3 lists, with one being
        # (stranded) PROseq GRanges, another being (unstranded) ATAC-seq data,
        # and the 3rd being (unstranded) ChIP-seq data. Each sublist will be independently
        # groupscaled based on the min/max of that data within the region.
        # List MUST be a named list, construct by calling:
        # list("PROseq" == PROseq.list,
        #      "ATACseq" == ATACseq.list, etc)
        # Parameters:
        # regions.gr: A GRanges object with a single gene. must have start, stop
        #             and strand, as well as one extra column called "tx_name"
        # dataset.grlist: A list of lists of GRanges objects. See above.
        # binsize: size of bins for applying bins_FUN over to aggregate signal.
        #           In general, less bins means less resolution but faster execution.
        # bins_FUN: function to be applied over bins (defaule = sum)
        # pad_left: integer of number of bp to add to right of gene
        # pad_right: integer of number of bp to add to left of gene
        # scale_bar_size: size of scale bar to include in bp
        # color_plus: value to use for coloring plus stranded signal
        # color_minus: value to use for coloring minus stranded signal
        # color_unstranded: value to use for coloring unstranded signal
        # tx_name: name of the mcol of regions.gr that contains the name
        #           of the transcript for plotting
        # labs: A character vector of labels for the signal tracks.
        #       Default is to use the names passed in through the lists.
        #       If you override, the length of the vector MUST match the total
        #       number of samples. The order of labels is the order of datasets
        #       in your list of signal tracks.
        # .expand_ranges: A named list of length = # data types, indicating
        #       whether each data type is in single width bins or run-compressed.
        #       Names of this list must be exact matches to names of your list of signal.
        #       If you aren't sure what this is, leave the default, which is
        #       to expand all data. If you know what you are doing, specifying
        #       FALSE for each datatype that is already in single base intervals
        #       speeds up getting counts ~2-fold.
       
        require(BRGenomics)
        require(ggplot2)
        require(facetscales)
        require(gtable)
        require(grid)
        require(readr)
        require(purrr)
        require(dplyr)
        
        # Main function for getting signal dataframes
        .get_signal_for_browser.jj <- 
            function(signal_name){
            # This function takes a single region as a GRanges object, and
            # A list of signal GRanges. Returns a dataframe for plotting
            # 'browser shots' with ggplot, with columns for count, coord,
            # sample (for faceting), and strand
            
            # Determining if all data in dataset.grlist is stranded or unstranded
            # Getting vector of unique values in the strand 
            # field across all GR objects in list
                
            .get_strand <- function(gr){
                # Input: GRanges object
                # Output: Chr vector of all unique values in the strand column
                return(unique(strand(gr)))
            }
            
            strand.chr <-
                as.character(unique(mclapply(dataset.grlist[[signal_name]], .get_strand))[[1]])
            
            # initializing values
            stranded = FALSE
            unstranded = FALSE
            
            # If stranded, setting stranded = T
            if("+" %in% strand.chr & "-" %in% strand.chr & !("*" %in% strand.chr)) {
                stranded = TRUE
            }
            
            # If unstranded, setting unstranded = T
            if(!("+" %in% strand.chr) & !("-" %in% strand.chr) & "*" %in% strand.chr) {
                unstranded = TRUE
            }
            
            # Throwing stop error if stranded or unstranded isn't set
            if (!(stranded | unstranded)) {
                stop(
                    "Your GRanges data isn't consistantly stranded/unstranded. It is okay to have a list of lists, but each sublist must have only the same type of data"
                )
            }
            
            # Function for getting list of signal matrices
            # TODO implement with future_map or mclapply for speed
            .get_signal_list <- function(regions.gr) {
                return(
                    getCountsByPositions(
                        dataset.grlist[[signal_name]],
                        regions.gr,
                        binsize = binsize,
                        expand_ranges = .expand_ranges[[signal_name]],
                        FUN = bin_FUN
                ))
            }
            
            # Function for generating binned signal dataframes coordinates from a 
            # vector of scores
            .make_signal_df <- function(sample_name, signal.lst, regions.gr){
                
                # If regions passed in are on the minus strand, the order is reveresed
                if(as.character(strand(regions.gr)) == "-"){
                    signal.lst[[sample_name]] <- rev(signal.lst[[sample_name]])
                }
                if(as.character(strand(regions.gr)) == "+"){
                    signal.lst[[sample_name]] <- t(signal.lst[[sample_name]])
                }
            
                return(
                    data.frame(
                        "count" = signal.lst[[sample_name]],
                        "coord" = seq(start(regions.gr), 
                                      end(regions.gr), 
                                      length.out = length(signal.lst[[sample_name]])),
                        "sample" = sample_name,
                        "strand" = ifelse(
                            unstranded,
                            "*",
                            as.character(strand(regions.gr))
                    )
                ))
            }
            
            # Generating dataframe of signal
            if(stranded){
                # generate plus and minus copies of regions.gr
                regions_plus.gr <- regions.gr
                regions_minus.gr <- regions.gr
                strand(regions_plus.gr) <- "+"
                strand(regions_minus.gr) <- "-"

                # Getting lists of signal
                signal_plus.lst <- .get_signal_list(regions_plus.gr)
                signal_minus.lst <- .get_signal_list(regions_minus.gr)


                # Getting df
                signal.df <- rbind(
                    map_dfr(
                        .x = names(dataset.grlist[[signal_name]]),
                        .f = .make_signal_df,
                        signal.lst = signal_plus.lst,
                        regions.gr = regions_plus.gr),
                    map_dfr(
                        .x = names(dataset.grlist[[signal_name]]),
                        .f = .make_signal_df,
                        signal.lst = signal_minus.lst,
                        regions.gr = regions_minus.gr))
            }
            
            if(unstranded){
                signal.lst <- .get_signal_list(regions.gr)
                
                # Getting df
                signal.df <-
                    map_dfr(
                        .x = names(dataset.grlist[[signal_name]]),
                        .f = .make_signal_df,
                        signal.lst = signal.lst,
                        regions.gr = regions.gr)
            }

            # Negating counts of minus strand signal
            if(stranded){
                signal.df$count[which(signal.df$strand == "-")] <-  
                    signal.df$count[which(signal.df$strand == "-")] * -1
            }
            
            return(signal.df)
        }
        
        # Custom theme for ggplot
        .ggtheme_browsershot.jj <- function() {
            theme_bw(base_size=10, base_family="Helvetica") %+replace%
                theme(
                    axis.text = element_text(size = 8), 
                    axis.ticks = element_line(colour = "black"), 
                    axis.line.y = element_line(color = "black"),
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.text.x = element_blank(),
                    legend.key = element_blank(), 
                    panel.background = element_blank(), 
                    panel.border = element_blank(), 
                    panel.grid.major.x = element_blank(), 
                    panel.grid.minor.x = element_blank(), 
                    panel.grid.major.y = element_blank(), 
                    panel.grid.minor.y = element_blank(), 
                    strip.background = element_blank(),
                    strip.text = element_text(size=10),
                    strip.text.y = element_text(angle = 360),
                    panel.spacing = unit(3, "mm")
                )
        }
        
   
        # Check that if labs were provided, there are the correct number
        if(!(length(.labs) == length(as.character(unlist(
            lapply(dataset.grlist, names)
        ))))) {
            stop("The vector of labels you provided is not the correct length")
        }
        
        # Check that length of expand_ranges == number of data types
        if(!(any(names(.expand_ranges) == names(dataset.grlist)))){
            stop("The logical list you passed to .expand_values does not have the same names as your dataset.grlist")
        }
        
        # Check that length of regions.gr == 1
        if(!(length(regions.gr) == 1)){
            stop("Your regions.gr does not include exactly one gene. Plotting multiple genes not currently supported.")
        }
        
        # Expanding regions.gr by pad but keeping original for annotations later
        annot.gr <- regions.gr
        start(regions.gr) <- start(regions.gr) - pad_left
        end(regions.gr) <- end(regions.gr) + pad_right
        
        # Getting center of both intervals for plotting annotation layer
        center <- start(regions.gr) + ((end(regions.gr) - start(regions.gr)) / 2 )
        center.annot <- start(annot.gr) + ((end(annot.gr) - start(annot.gr)) / 2 )
        
        # Getting signal matrix for each element of dataset.grlist
        signal.lst <- mclapply(X = names(dataset.grlist), FUN = .get_signal_for_browser.jj)

        
        # Functions for programatically determining yscales
        # with groupscale per list element of dataset.grlist
        .get_minmax <- function(df,
                                field = "count",
                                enforce0 = TRUE) {
            # takes a df as input and returns a vector of c(min, max)
            # of the range of the data, ceilinged at the top and floored
            # at the bottom (ie -5.1 to 10.1 becomes c(-6, 11)).
            # Useful for programatically defining the scale breaks of a plot
            # if enforce0 == TRUE (default), if the data is entirely positive, the
            # min is fixed at 0 and if it's all negative the max is fixed at 0
            # ie 5.4-10.1 becomes (0, 11) and -12.3 to -4.5 becomes (-13,0)
            ymin <- floor(min(df[field]))
            ymax <- ceiling(max(df[field]))
            if (enforce0) {
                if (ymin > 0) {
                    ymin <- 0
                }
                if (ymax < 0) {
                    ymax <- 0
                }
            }
            if(
                (ymin < 0 & ymax > 0) & (0-ymin > ymax*0.1) & !((0-ymin)*0.1 > ymax)
                ){
                return(c(ymin, 0, ymax))
            } else {
                return(c(ymin, ymax))
            }
        }
        
        .get_yscale <- function(df, field = "sample") {
            # Takes a list of signal dfs and generates a list of
            # y scales for each unique value in field (facetting variable)
            vals <- .get_minmax(df, field ="count")
            scale <- scale_y_continuous(
                limits = c(vals[1], vals[length(vals)]),
                breaks = vals,
                labels = vals,
                expand = c(0, 0)
            )
            scale.lst <- c()
            for (i in unique(df[[field]])) {
                scale.lst[[i]] = scale
            }
            return(scale.lst)
        }
        
        
        
        # Getting scales for each data type
        scales_y <- flatten(lapply(signal.lst, .get_yscale))
            
        
        # Adding y_scale for genes panel 
        # (numbers are arbitrary and just for anchoring objects)
        scales_y$genes <- scale_y_continuous(limits = c(0, 5),
                                             breaks = NULL,
                                             labels = NA)
        
        # unlisting signal
        signal.df <- do.call("rbind", signal.lst)
        
        # Converting sample column to chr
        signal.df$sample <- as.character(signal.df$sample)
        
        # Adding single line to create gene annotation level
        signal.df <- rbind(
            signal.df,
            c(NA, min(signal.df$coord), "genes", "+", "genes"))
        
        # Factoring signal (in order it was in lists) with genes at top level
        signal.df$sample <- factor(
            signal.df$sample, 
            levels = c("genes", 
                as.character(unlist(lapply(dataset.grlist, names)))))
        
        # Setting strand to fwd, rev, or all
        signal.df <-
            mutate(signal.df, strand = ifelse(
                strand == "+", "fwd", ifelse(strand == "-", "rev", "all")))
        
        # Factoring strand
        signal.df$strand <-
            factor(signal.df$strand, levels = c("fwd", "rev", "all"))
        
        # Converting count and coord to numeric
        signal.df$count <- as.numeric(signal.df$count)
        signal.df$coord <- as.numeric(signal.df$coord)
        
        # Creating df of features
        features.df <- 
            data.frame(
                "type" = c("gene", "scalebar", "chr"),
                "xmin" = c(start(annot.gr), center - 0.5*scale_bar_size, center),
                "xmax" = c(end(annot.gr), center + 0.5*scale_bar_size, center),
                "ymin" = c(0, 2.8, 4.5),
                "ymax" = c(1.8, 2.8, 4.5),
                "sample" = c("genes", "genes", "genes"),
                "strand" = rep(as.character(strand(regions.gr)), 3),
                "label" = NA,
                "label_center" = rep(center, 3)
            )
        
        # Gene label
        features.df$label[1] <-
            ifelse(as.character(strand(regions.gr)) == "-",
                   paste0("<<< ", regions.gr$tx_name),
                   paste0(regions.gr$tx_name, " >>>"))
        
        # scalebar label
        features.df$label[2] <- 
            paste0(as.character(scale_bar_size /1000), "kb")
        
        # chr label
        features.df$label[3] <- 
            paste0(as.character(seqnames(regions.gr)),
                   ":",
                   as.character(start(regions.gr)),
                   "-",
                   as.character(end(regions.gr)))
        
        # Set 0s to NA (othersiwe gglot draws lines at 0
        # that look like low signal where there actually is nothing)
        signal.df$count[which(signal.df$count == 0)] <- NA
        
        
        # Creating labeller (inserting "" so that the annotation track isn't labelled)
        .labs <- c("", .labs)
        names(.labs) <- levels(signal.df$sample)
        
    
        
        p <- ggplot(signal.df, aes(x = coord, y = count)) +
            geom_col(show.legend = F, aes(fill = strand, color = strand), size = 0.5)+
            .ggtheme_browsershot.jj() +
            facet_grid_sc(rows = vars(sample),
                         scales = list(y = scales_y),
                         labeller = labeller(sample = .labs))+
            scale_color_manual(values = c(
                "fwd" = color_plus,
                "rev" = color_minus,
                "all" = color_unstranded
            )) +
            scale_fill_manual(values = c(
                "fwd" = color_plus,
                "rev" = color_minus,
                "all" = color_unstranded
            )) +
            scale_x_continuous(expand = c(0,0))+
            geom_rect(
                data = features.df[1,],
                inherit.aes = F,
                aes(
                    xmin = xmin,
                    xmax = xmax,
                    ymin = ymin,
                    ymax = ymax
                )
            ) +
            geom_text(
                data = features.df[1,],
                inherit.aes = F,
                aes(x = label_center, y = 0.9, label = label),
                color = "white",
                fontface = 3,
                size = 3
            )+
            geom_segment(
                data = features.df[2,],
                inherit.aes = F,
                aes(
                    x = xmin,
                    xend = xmax,
                    y = ymin,
                    yend = ymax)
                )+
            geom_text(
                data = features.df[2,],
                inherit.aes = F,
                aes(x = (xmin - 0.25*scale_bar_size), y = ymin, label = label),
                color = "black",
                size = 3
            )+
            geom_text(
                data = features.df[3,],
                inherit.aes = F,
                aes(x = label_center, y = ymin, label = label),
                color = "black",
                size = 3
            )+
            xlab(NULL)+
            ylab(NULL)

        # getting gtable from plot
        p_tab <- ggplotGrob(p)
        
        # Removing left axis 1 (the one in the gene annotation track)
        .filter_gtable <- function(gtable, drop.chr){
            # Takes a gtable object and removes grobs in drop.
            # dropcan be chr or a vector of chr.
            hits <- !(gtable$layout$name %in% drop.chr) # Get grobs not in drop
            gtable$layout <- gtable$layout[hits, , drop = FALSE]
            gtable$grobs <- gtable$grobs[hits]
            gtable <- gtable_trim(gtable)
        }
        
        p_tab <- .filter_gtable(p_tab, "axis-l-1")
        return(p_tab)
    }
