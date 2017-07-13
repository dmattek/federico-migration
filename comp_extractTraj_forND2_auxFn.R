require(data.table)
require(ggplot2)
require(ggrepel) # for non-overlapping labels

# definition of custom colour palette
rhg_cols <- c("#771C19","#AA3929","#E25033","#F27314","#F8A31B",
              "#E2C59F","#B6C5CC","#8E9CA3","#556670","#000000")

md_cols <- c("#FFFFFF",
             "#F8A31B","#F27314",  "#E25033",  "#AA3929",
             "#FFFFCC", "#C2E699", "#78C679", "#238443")


myCheckDigits <- function(x){ grepl('^[-]?[0-9]+[.]?[0-9]*$' , x) }

myCheckLogical <- function(x){ grepl('^TRUE$|^FALSE$' , x) }

myAtan = function(x, y) {
  z = 0
  
  if (is.na(x) || is.na(y))
    z = NA_real_
  else {
    if (x > 0 && y > 0) 
      z = atan(x/y)
    
    if (x > 0 && y == 0)
      z = pi * 0.5
    
    if (x > 0 && y < 0)
      z = atan(abs(y)/x) + pi * 0.5
    
    if (x == 0 && y < 0)
      z = pi
    
    if (x < 0 && y < 0)
      z = atan(x/y) + pi
    
    if (x < 0 && y == 0)
      z = pi * 1.5
    
    if (x < 0 && y > 0)
      z = atan(y/abs(x)) + pi * 1.5
    
  }
  
  return(z)
}

myRotate = function(x, y, th) {
  list(
    x * cos(th) + y * sin(th),
    -x * sin(th) + y * cos(th)
  )
}

myConvertStringListToTypes <- function(in.l){ 
  # convert strings with digits to numeric
  # uses logical indexing: http://stackoverflow.com/questions/42207235/replace-list-elements-by-name-with-another-list
  loc.l = myCheckDigits(in.l)
  in.l[loc.l] = lapply(in.l[loc.l], as.numeric)
  
  # convert strings with TRUE/FALSE to logical
  loc.l = myCheckLogical(in.l)
  in.l[loc.l] = lapply(in.l[loc.l], as.logical)
  
  return(in.l) 
}

myTrajExtr = function(in.dt,
                      in.max.break = 1, 
                      in.met.series = 'Metadata_Series', 
                      in.met.t = 'Metadata_T', 
                      in.met.tracklabel = 'TrackObjects_Label',
                      in.aggr = FALSE,
                      in.aggr.cols = NULL) {
  
  loc.dt = copy(in.dt)
  
  # TrackObjects assigns the same label for different cells from IdentifyPrimaryObjects
  # The following aggregation makes sure there's a unique TrackObjects_Label
  # for every Site at every time point: it takes the mean intensity of duplicated labels.
  # Make sure it makes sense!!!
  # Roughly 10% of objects affected
  
  if (in.aggr) {
    loc.dt = loc.dt[, lapply(.SD, mean),  
                    by = c(in.met.series, in.met.t, in.met.tracklabel), .SDcols = in.aggr.cols]
  }
  
  
  # cells from different sites have the same TrackObjects_Label
  # make it unique accross the experiment and pad numbers with zeros
  loc.dt[, TrackObjects_Label_uni := paste(sprintf("%02d", get(in.met.series)), sprintf("%04d", get(in.met.tracklabel)), sep = "_")]
  
  # set the length of trajectories to be selected
  loc.n.tp = max(loc.dt[, in.met.t, with = FALSE]) - min(loc.dt[, in.met.t, with = FALSE])
  
  # Select only cells with at least 1st and last timepoint measured
  # Metadata_T starts at 0, thus cells with T=0 & T=last measured 
  # have max(Metadata_T) - min(Metadata_T) = n.tp
  
  loc.dt.tmp1 = loc.dt[, .(Metadata_T.span = max(get(in.met.t)) - min(get(in.met.t))), 
                       by = TrackObjects_Label_uni][Metadata_T.span == loc.n.tp]
  
  
  # With cells selected above,
  # select tcourses with at most 1 break.
  ## Single-breaks usually result from a single frame with lost focus
  ## no interpolation
  # Calculation based on differences between consecutive Metadata_T.
  # Metadata_T increases by 1, thus a single-frame break would result in a difference of 2.
  
  loc.dt.tmp2 = loc.dt[TrackObjects_Label_uni %in% loc.dt.tmp1$TrackObjects_Label_uni]
  loc.dt.tmp2[, Metadata_T.diff := c(NA, diff(get(in.met.t))), by = TrackObjects_Label_uni]
  loc.dt.tmp2 = loc.dt.tmp2[Metadata_T.diff <= in.max.break + 1]
  
  # Selected trajectories with at most 1-frame break
  loc.out = loc.dt[TrackObjects_Label_uni %in% loc.dt.tmp2$TrackObjects_Label_uni]
  return(loc.out)
}


myTrajExtr2 = function(in.dt,
                       in.max.break = 1, 
                       in.met.series = 'Metadata_Series', 
                       in.met.t = 'Metadata_T', 
                       in.met.tracklabel = 'TrackObjects_Label') {
  loc.dt = copy(in.dt)
  
  loc.t.range = unique(loc.dt[, c(in.met.t), with = FALSE])
  
  loc.dt[, TrackObjects_Label_uni := paste(sprintf("%02d", get(in.met.series)), sprintf("%04d", get(in.met.tracklabel)), sep = "_")]
  loc.dt.tmp1 = loc.dt[, .(Ntpt = .N, 
                           T.start = first(Metadata_T), 
                           T.end = last(Metadata_T)), 
                       by = TrackObjects_Label_uni][Ntpt <=  nrow(loc.t.range) & 
                                                      T.start == min(loc.t.range) & 
                                                      T.end == max(loc.t.range)]
  
  
  # With cells selected above,
  # select tcourses with at most 1 break.
  ## Single-breaks usually result from a single frame with lost focus
  ## no interpolation
  # Calculation based on differences between consecutive Metadata_T.
  # Metadata_T increases by 1, thus a single-frame break would result in a difference of 2.
  
  loc.dt.tmp2 = loc.dt[TrackObjects_Label_uni %in% loc.dt.tmp1$TrackObjects_Label_uni]
  loc.dt.tmp2[, Metadata_T.diff := c(NA, diff(get(in.met.t))), by = TrackObjects_Label_uni]
  loc.dt.tmp2 = loc.dt.tmp2[Metadata_T.diff <= in.max.break + 1]
  
  # Selected trajectories with at most 1-frame break
  loc.out = loc.dt[TrackObjects_Label_uni %in% loc.dt.tmp2$TrackObjects_Label_uni]
  return(loc.out)
}


myTrajExtr3 = function(in.dt,
                       in.max.break = 1, 
                       in.met.t = 'Metadata_T', 
                       in.met.tracklabel = 'TrackObjects_Label') {
  loc.dt = copy(in.dt)
  
  loc.t.range = unique(loc.dt[, c(in.met.t), with = FALSE])
  
  loc.dt.tmp1 = loc.dt[, .(Ntpt = .N, 
                           T.start = first(Metadata_T), 
                           T.end = last(Metadata_T)), 
                       by = in.met.tracklabel][Ntpt <=  nrow(loc.t.range) & 
                                                      T.start == min(loc.t.range) & 
                                                      T.end == max(loc.t.range)]
  
  
  # With cells selected above,
  # select tcourses with at most 1 break.
  ## Single-breaks usually result from a single frame with lost focus
  ## no interpolation
  # Calculation based on differences between consecutive Metadata_T.
  # Metadata_T increases by 1, thus a single-frame break would result in a difference of 2.
  
  loc.dt.tmp2 = loc.dt[get(in.met.tracklabel) %in% loc.dt.tmp1[, get(in.met.tracklabel)]]
  loc.dt.tmp2[, Metadata_T.diff := c(NA, diff(get(in.met.t))), by = in.met.tracklabel]
  loc.dt.tmp2 = loc.dt.tmp2[Metadata_T.diff <= in.max.break + 1]
  
  # Selected trajectories with at most 1-frame break
  loc.out = loc.dt[ get(in.met.tracklabel) %in% loc.dt.tmp2[, get(in.met.tracklabel)]]
  return(loc.out)
}


# Returns original dt with an additional column with normalized quantity
# the column to be normalised is given by 'in.meas.col'
# Normalisation is based on part of the trajectory;
# this is defined by in.rt.in and max, and the column with time in.rt.col.
# Additional parameters:
# in.by.cols - character vector with 'by' columns to calculate normalisation per group
#              if NULL, no grouping is done
# in.robust - whether robust measures should be used (median instead of mean, mad instead of sd)
# in.type - type of normalization: z.score or mean (fold change w.r.t. mean)

myNorm = function(in.dt,
                  in.meas.col,
                  in.rt.col = 'RealTime',
                  in.rt.min = 10,
                  in.rt.max = 20,
                  in.by.cols = NULL,
                  in.robust = TRUE,
                  in.type = 'z.score') {
  loc.dt <- copy(in.dt) # copy so as not to alter original dt object w intermediate assignments
  
  if (is.null(in.by.cols)) {
    if (in.robust)
      loc.dt.pre.aggr = loc.dt[get(in.rt.col) > in.rt.min & get(in.rt.col) < in.rt.max, .(meas.md = median(get(in.meas.col)),
                                                                                          meas.mad = mad(get(in.meas.col)))]
    else
      loc.dt.pre.aggr = loc.dt[get(in.rt.col) > in.rt.min & get(in.rt.col) < in.rt.max, .(meas.md = mean(get(in.meas.col)),
                                                                                          meas.mad = sd(get(in.meas.col)))]
    
    loc.dt = cbind(loc.dt, loc.dt.pre.aggr)
  }  else {
    if (in.robust)
      loc.dt.pre.aggr = loc.dt[get(in.rt.col) > in.rt.min & get(in.rt.col) < in.rt.max, .(meas.md = median(get(in.meas.col)),
                                                                                          meas.mad = mad(get(in.meas.col))), by = in.by.cols]
    else
      loc.dt.pre.aggr = loc.dt[get(in.rt.col) > in.rt.min & get(in.rt.col) < in.rt.max, .(meas.md = mean(get(in.meas.col)),
                                                                                          meas.mad = sd(get(in.meas.col))), by = in.by.cols]
    
    loc.dt = merge(loc.dt, loc.dt.pre.aggr, by = in.by.cols)
  }
  
  
  if (in.type == 'z.score') {
    loc.dt[, meas.norm := (get(in.meas.col) - meas.md) / meas.mad]
  } else {
    loc.dt[, meas.norm := (get(in.meas.col) / meas.md)]
  }
  
  loc.dt[, c('meas.md', 'meas.mad') := NULL]
  return(loc.dt)
}


# add a rescaled column by shifting a selected column to a baseline (in.y.min - in.y.h) and setting height to in in.y.h
# the column to shift is in.meas.col
# the new column has the same name as the original but with a suffix in.meas.resc.col
myResc = function(in.dt,
                  in.meas.col,
                  in.meas.resc.col = '.resc',
                  in.y.min,
                  in.y.h) {
  
  loc.dt = copy(in.dt)
  loc.dt[, paste0(in.meas.col, in.meas.resc.col) := (in.y.min - in.y.h) +
           (in.y.min - (in.y.min - in.y.h)) * (get(in.meas.col) - min(get(in.meas.col))) /
           (max(get(in.meas.col)) - min(get(in.meas.col)))]
  
  return(loc.dt)
}

myWilcoxTest = function(in.data.x, 
                        in.data.cat, 
                        in.method = 'fdr', 
                        in.cut = c(0, 0.0001, 0.001, 0.01, 0.05, 1)) {
  test.res = pairwise.wilcox.test(in.data.x, in.data.cat, p.adjust.method = in.method, paired = FALSE, alternative = 'two.sided')
  test.res.pval = as.data.table(melt(test.res$p.value))
  test.res.pval[, pval.levels := cut(value, in.cut, right = TRUE, include.lowest = TRUE) ]
  
  return(test.res.pval)
}



myGgplotTraj = function(dt.arg,
                        x.arg,
                        y.arg,
                        group.arg,
                        facet.arg = NULL,
                        summary.arg = FALSE,
                        facet.ncol.arg = 2,
                        line.col.arg = NULL,
                        xlab.arg = "Time",
                        ylab.arg = "Fl. int.",
                        plotlab.arg = "",
                        dt.stim.arg = NULL,
                        stim.x.arg,
                        stim.y.arg,
                        maxrt.arg = 60,
                        xaxisbreaks.arg = 10,
                        xlim.arg = NULL,
                        ylim.arg = NULL) {
  
  p.tmp = ggplot(dt.arg,
                 aes_string(x = x.arg,
                            y = y.arg,
                            group = group.arg))
  
  if (is.null(line.col.arg))
    p.tmp = p.tmp + geom_line(alpha = 0.25, size = 0.25)
  else
    p.tmp = p.tmp + geom_line(aes_string(colour = line.col.arg), alpha = 0.5, size = 0.5)
  
  
  if (summary.arg)
    p.tmp = p.tmp + 
      stat_summary(aes_string(y = y.arg, group = 1), 
                   geom="ribbon", 
                   fun.data=mean_cl_boot, 
                   colour='red', alpha=0.5,
                   group = 1) +
      stat_summary(
        aes_string(y = y.arg, group = 1),
        fun.y = mean,
        colour = 'red',
        linetype = 'solid',
        size = 1,
        geom = "line",
        group = 1
      )
  
  
  
  if (!is.null(facet.arg))
    p.tmp = p.tmp + 
      facet_wrap(as.formula(paste("~", facet.arg)),
                 ncol = facet.ncol.arg,
                 scales = "free_x", 
                 drop = FALSE)
  
  if(!is.null(dt.stim.arg)) {
    p.tmp = p.tmp + geom_line(data = dt.stim.arg,
                              aes_string(x = stim.x.arg, y = stim.y.arg),
                              colour = 'blue',
                              size = 1,
                              group = 1) 
  }
  
  if (!is.null(ylim.arg))
    p.tmp = p.tmp +
      coord_cartesian(ylim = ylim.arg)
  
  if (!is.null(xlim.arg))
    p.tmp = p.tmp +
      coord_cartesian(xlim = xlim.arg)
  
  p.tmp = p.tmp + 
    scale_x_continuous(breaks = seq(0, maxrt.arg, xaxisbreaks.arg)) +
    xlab(paste0(xlab.arg, "\n")) +
    ylab(paste0("\n", ylab.arg)) +
    ggtitle(plotlab.arg) +
    theme_bw(base_size = 18, base_family = "Helvetica") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black", size = 0.25),
      axis.line.y = element_line(color = "black", size = 0.25),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      strip.text.x = element_text(size = 14, face = "bold"),
      strip.text.y = element_text(size = 14, face = "bold"),
      strip.background = element_blank(),
      legend.key = element_blank(),
      legend.key.height = unit(1, "lines"),
      legend.key.width = unit(2, "lines"),
      legend.position = "top"
    )
  
  return(p.tmp)
}

myGgplotTrajRibbon = function(dt.arg,
                              x.arg,
                              y.arg,
                              group.arg, 
                              xlab.arg = "Time",
                              ylab.arg = "Fl. int.",
                              plotlab.arg = "") {
  p.tmp = ggplot(dt.arg, aes_string(x = x.arg, group = group.arg)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = 'grey70', alpha = 0.5) +
    geom_line(aes_string(y = y.arg, colour = group.arg)) +
    scale_color_manual(values = rhg_cols, name = "") +
    xlab(paste0(xlab.arg, "\n")) +
    ylab(paste0("\n", ylab.arg)) +
    ggtitle(plotlab.arg) +
    theme_bw(base_size = 18, base_family = "Helvetica") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black", size = 0.25),
      axis.line.y = element_line(color = "black", size = 0.25),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      strip.text.x = element_text(size = 14, face = "bold"),
      strip.text.y = element_text(size = 14, face = "bold"),
      strip.background = element_blank(),
      legend.key = element_blank(),
      legend.key.height = unit(1, "lines"),
      legend.key.width = unit(2, "lines"),
      legend.position = "top"
    ) 
  return(p.tmp)
}

myGgplotXYtraj = function(dt.arg,
                          x.arg,
                          y.arg,
                          group.arg,
                          x.rev.arg = FALSE,
                          y.rev.arg = FALSE,
                          facet.arg = NULL,
                          facet.ncol.arg = 2,
                          xlab.arg = "Position X (pixels)",
                          ylab.arg = "Position Y (pixels)",
                          plotlab.arg = "",
                          path.size.arg = 1,
                          path.alpha.arg = 0.5,
                          path.col.arg = NULL,
                          path.col.pal.arg = "RdYlBu",
                          point.arg = NULL,
                          point.size.arg = 5,
                          point.col.arg = 'position') {
  
  p.tmp = ggplot(dt.arg, aes_string(x = x.arg, y = y.arg))
  
  if (is.null(path.col.arg)) {
    p.tmp = p.tmp +
      geom_path(aes_string(group = group.arg), 
                size = path.size.arg, 
                alpha = path.alpha.arg, 
                lineend = "round")
    
  } else {
    p.tmp = p.tmp +
      geom_path(aes_string(group = group.arg, colour = path.col.arg), 
                size = path.size.arg, 
                alpha = path.alpha.arg, 
                lineend = "round") +
      scale_color_distiller(palette = "RdYlBu", 
                            name = 'Direction (rad)',
                            limit = pi*c(0, 2), 
                            breaks = pi*seq(0, 2, .5), 
                            labels = c('0', 'Pi/2', 'Pi', '3/2 Pi', '2 Pi'))
  }
  
  if (x.rev.arg)
    p.tmp = p.tmp +
      scale_x_reverse()
  
  if (y.rev.arg)
    p.tmp = p.tmp +
      scale_y_reverse()
  
  if (!is.null(point.arg))
    p.tmp = p.tmp +
      geom_point(data = point.arg, 
                 aes_string(colour = point.col.arg), 
                 size = point.size.arg) +
      scale_colour_manual(values = c('green', 'red'))
  
  if (!is.null(facet.arg))
    p.tmp = p.tmp +
      facet_wrap(
        as.formula(paste("~", facet.arg)),
        ncol = facet.ncol.arg,
        scales = "fixed",
        drop = FALSE
      )
  
  p.tmp = p.tmp +
    xlab(paste0(xlab.arg, "\n")) +
    ylab(paste0("\n", ylab.arg)) +
    ggtitle(plotlab.arg) +
    theme_bw(base_size = 18, base_family = "Helvetica") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black", size = 0.25),
      axis.line.y = element_line(color = "black", size = 0.25),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      strip.text.x = element_text(size = 14, face = "bold"),
      strip.text.y = element_text(size = 14, face = "bold"),
      strip.background = element_blank(),
      legend.key = element_blank(),
      legend.key.height = unit(1, "lines"),
      legend.key.width = unit(2, "lines"),
      legend.position = "top"
    )
  
  return(p.tmp)
}


myGgplotDens = function(dt.arg,
                        x.arg,
                        group.arg,
                        xlab.arg = "Time",
                        ylab.arg = "Fl. int.",
                        plotlab.arg = "",
                        vline.dt.arg = NULL) {
  p.tmp = ggplot(dt.arg, aes_string(x = x.arg), group = group.arg) +
    geom_density(aes_string(y = "..density..", colour = group.arg), size = .5)
  
  if (!is.null(vline.dt.arg))
    p.tmp = p.tmp + 
      geom_vline(data = vline.dt.arg, 
                 aes_string(xintercept = 'V1', 
                            group = group.arg, 
                            colour = group.arg), 
                 linetype="dashed") +
      geom_text_repel(data = vline.dt.arg, aes( V1, 0, label = sprintf("%.2f", V1)), size = 3)
  
  p.tmp = p.tmp + 
    scale_color_discrete(name = "") +
    xlab(paste0(xlab.arg, "\n")) +
    ylab(paste0("\n", ylab.arg)) +
    ggtitle(plotlab.arg) +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black", size = 0.25),
      axis.line.y = element_line(color = "black", size = 0.25),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      strip.text.x = element_text(size = 12, face = "bold"),
      strip.text.y = element_text(size = 12, face = "bold"),
      strip.background = element_blank(),
      legend.key = element_blank(),
      legend.key.height = unit(1, "lines"),
      legend.key.width = unit(2, "lines"),
      legend.position = "top"
    )
  return(p.tmp)
}


myGgplotBoxPlot = function(dt.arg,
                           x.arg,
                           y.arg,
                           group.arg,
                           xlab.arg = "Time",
                           ylab.arg = "Fl. int.",
                           plotlab.arg = "",
                           y.lim.perc = c(0.05, 0.95),
                           boxplot.notch = TRUE,
                           dotplot.view = TRUE,
                           Dots_Dim = 0.01) {
  
  p.tmp = ggplot(dt.arg, aes_string(x.arg, y.arg)) +
    geom_boxplot(outlier.shape = NA, notch = boxplot.notch)
  
  if(dotplot.view)
    p.tmp = p.tmp + geom_dotplot(binwidth = Dots_Dim , binaxis = "y", stackdir = "center", 
                 fill = "black", color = "black")
  
    p.tmp = p.tmp + coord_cartesian(ylim = quantile(dt.arg[, get(y.arg)], y.lim.perc, na.rm = TRUE)) +
    scale_color_discrete(name = "") +
    xlab(paste0(xlab.arg, "\n")) +
    ylab(paste0("\n", ylab.arg)) +
    ggtitle(plotlab.arg) +
    theme_bw(base_size = 10, base_family = "Helvetica") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black", size = 0.25),
      axis.line.y = element_line(color = "black", size = 0.25),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      strip.text.x = element_text(size = 12, face = "bold"),
      strip.text.y = element_text(size = 12, face = "bold"),
      strip.background = element_blank(),
      legend.key = element_blank(),
      legend.key.height = unit(1, "lines"),
      legend.key.width = unit(2, "lines"),
      legend.position = "top"
    )
  return(p.tmp)
}


myRasterPlot = function(in.data, 
                        in.x = 'Var1', 
                        in.y = 'Var2', 
                        in.fill = 'pval.levels', 
                        in.title.string = '',
                        in.subtitle.string = '',
                        in.legend.string = 'p-value:',
                        in.breaks = rev(c( '[0,0.0001]', '(0.0001,0.001]', '(0.001,0.01]', '(0.01,0.05]', '(0.05,1]','<NA>')),
                        in.labels = rev(c('<= 0.0001', '<= 0.001', "<= 0.01", '<= 0.05', '> 0.05',"")),
                        in.col.vals = rhg_cols[(c(3,4,5,7,8))]) {
  
  p.out = ggplot(test.res.pers, aes_string(x = in.x, y = in.y, fill = in.fill)) + 
    geom_raster() + 
    scale_fill_manual(name = in.legend.string,
                      na.value = 'white',
                      breaks = in.breaks, 
                      values = in.col.vals, 
                      label = in.labels, 
                      drop = FALSE) +
    xlab("") +
    ylab("") +  
    ggtitle(in.title.string, subtitle = in.subtitle.string) +
    theme_bw(base_size = 18, base_family = "Helvetica") +
    theme(panel.grid.minor=element_blank(), 
          panel.grid.major=element_blank(),
          panel.border=element_blank(),
          axis.line = element_line(colour="black"),
          axis.text.x = element_text(size=10, angle = 45, hjust = 1),
          strip.text.x = element_text(size=14, face="bold"),
          strip.text.y = element_text(size=14, face="bold"),
          strip.background = element_blank(),
          legend.key = element_blank(),
          legend.key.height = unit(1, "lines"),
          legend.key.width = unit(2, "lines"),
          legend.position = "right")
  
  return(p.out)
}
