####
## Analyzes Image.csv & Nuclei.csv
## Extracts single-cell trajectories

# Execute in a directory where cp.out is located
# This version handles multi-well experiments
#
# Image.csv should include:
#   ImageNumber
#   Metadata_Series
#   Metadata_T	
#   Metadata_Well	
#   Metadata_Experiment
#
# Nuclei.csv should include:
#   ImageNumber
#   TrackObjects_Label


#source(pipe(paste("wget -O -", "https://www.dropbox.com/s/qpk9wn9mfoujvog/comp_extractTraj_forND2_auxFn.R?dl=0")))

# loads required libraries
require(data.table)
require(dplyr)
require(ggplot2)
require(xlsx)
require(gridExtra)
require(Hmisc) # for smean.cl.boot
require(scales) # for rescale in gradient plotting
require(RCurl)


auxScript = getURL("https://raw.githubusercontent.com/dmattek/federico-migration/master/comp_extractTraj_forND2_auxFn.R", ssl.verifypeer = FALSE)
eval(parse(text = auxScript))


# definition of custom colour palette
rhg_cols <- c("#771C19","#AA3929","#E25033","#F27314","#F8A31B",
              "#E2C59F","#B6C5CC","#8E9CA3","#556670","#000000")

md_cols <- c("#FFFFFF",
             "#F8A31B","#F27314",  "#E25033",  "#AA3929",
             "#FFFFCC", "#C2E699", "#78C679", "#238443")


####
## Main parameters

# names of two files with paramaters for analysis and experiment
# these files should be located one folder up of cp.out from which the script is executed
s.fname.par = 'plotFormat.xlsx'

s.fname.exp = list.files(
  path = '../.',
  pattern = '*escription.xlsx',
  recursive = FALSE,
  full.names = TRUE
)

b.error = FALSE
if (length(s.fname.exp) > 1) {
  print("Error: mre than one file with experimental description found!")
  b.error = TRUE
}
stopifnot(!b.error)

####
## Read parameters from the file

df.par = read.xlsx(
  paste0('../', s.fname.par),
  sheetIndex = 1,
  header = FALSE,
  as.data.frame = TRUE,
  colClasses = rep("character", 2),
  colIndex = 1:2,
  stringsAsFactors = FALSE
)

# convert data frame with parameters to a list 
l.par = split(df.par[, 2], df.par[, 1])

# convert strings with digits to numeric and strings with TRUE/FALSE to logical
l.par = myConvertStringListToTypes(l.par)


#####
## Read experiment description
dt.exp = as.data.table(read.xlsx(s.fname.exp, 
                                 sheetIndex = 1, 
                                 header = TRUE, 
                                 startRow = 4, colIndex = if (l.par$exp.wound) 1:3 else 1:2,
                                 as.data.frame = TRUE, 
                                 stringsAsFactors=FALSE))

# sometimes an NA column appears at the end; remove
dt.exp = dt.exp[, names(dt.exp)[!(names(dt.exp) %like% 'NA')], with = FALSE]

# sometimes an NA row appears at the end; remove
dt.exp = dt.exp[dt.exp[, !Reduce(`&`, lapply(.SD, is.na))]]



####
## Read data

# dt.im file contains ImageNumber column and all metadata fields
# dt.nuc contains ImageNumber, measurement columns and TrackLabel (no metadata)
# For trajectory selection, a unique cell ID is created that consists of
# MetadataWell (if exists), MetadataSeries, and TrackObjectLabel.
# This identifies the cell in the entire experiment.

## Load Nuclei files using custom file reading function
dt.nuc = fread(paste(l.par$dir.out, l.par$files.nuc, sep = '/'))

# load Image files
dt.im = fread(paste(l.par$dir.out, l.par$files.img, sep = '/'))



## assignment of columns:
## s.met.site
## s.met.time
## s.tracklabel
## s.meas


# deciding whether FOV are in Image_Metadata_Site or Image_Metadata_Series
# assign to Image_Metadata_Site or Image_Metadata_Series
b.error = FALSE
if (length(names(dt.im)[names(dt.im) %like% 'Site']) == 1)
  s.met.site = names(dt.im)[names(dt.im) %like% 'Site'] else 
    if (length(names(dt.im)[names(dt.im) %like% 'Series']) == 1)
      s.met.site = names(dt.im)[names(dt.im) %like% 'Series'] else {
        print("Error extracting Site/Series metadata!")
        b.error = TRUE
      }
stopifnot(!b.error)

# assign to Image_Metadata_T or Image_Metadata_Time
s.met.time = names(dt.im)[names(dt.im) %like% 'Metadata_T']

# assign to TrackObjects_Label
s.met.tracklabel = names(dt.nuc)[names(dt.nuc) %like% 'TrackObjects_Label']

# experiment name column (used in plot titles)
s.exp.name.col = names(dt.im)[names(dt.im) %like% 'Exp']

# if more than 1 experiment in a file, use the 1st name only
s.exp.name = unique(dt.im[, get(s.exp.name.col)])
if(length(s.exp.name) > 1)
  s.exp.name = 'Experiment'

# assign image number column
s.met.imno = names(dt.im)[names(dt.im) %like% 'ImageNumber']


# If Metadata_Well exists, convert Site/Series to a unique combination of Well and Site/Series
s.met.well = names(dt.im)[names(dt.im) %like% 'Well']
b.well = FALSE
if (length(s.met.well) > 0) {
#  dt.im[, (s.met.site) := sprintf('%s_%02d', get(s.met.well), get(s.met.site))]

  b.well = TRUE
}

# Create a table with real time (in min) assigned to to Metadata_T
# Metadata_T starts at 0

dt.tass = data.table(Metadata_T = unique(dt.im[, get(s.met.time)]), 
                     RealTime.m = seq(0, max(unique(dt.im[, get(s.met.time)]))*l.par$exp.acquisition.freq, l.par$exp.acquisition.freq),
                     RealTime.h = seq(0, max(unique(dt.im[, get(s.met.time)]))*l.par$exp.acquisition.freq, l.par$exp.acquisition.freq) / 60)

# Trim to analysis duration
# Following merges will take into account this trimming
dt.tass = dt.tass[RealTime.m >= l.par$exp.skip & RealTime.m <= l.par$exp.analysis]


# Add Site/Series column to Nuclei (necessary for creating unique cell id in trajectory selection)
v.colnames = names(dt.im)
v.colnames = v.colnames[!(v.colnames %like% 'Experiment')]
dt.nuc = merge(dt.nuc, dt.im[, v.colnames, with = FALSE], by = s.met.imno)

# Create a column with unique cell ID
s.met.tracklabel.uni = paste0(s.met.tracklabel, '_uni')
if (b.well)
  dt.nuc[, (s.met.tracklabel.uni) := sprintf('%s_%02d_%04d', get(s.met.well), get(s.met.site), get(s.met.tracklabel))] else
  dt.nuc[, (s.met.tracklabel.uni) := sprintf('%02d_%04d', get(s.met.site), get(s.met.tracklabel))]


# add real time column
# this is only to trim the track to [exp.skip, exp.analysis] range
dt.nuc = merge(dt.nuc, dt.tass, by = s.met.time)

# obtain data table with cells with trajectories spanning entire experiment
if(b.well)
  setkeyv(dt.nuc, c(s.met.well, s.met.site, s.met.tracklabel.uni, s.met.time)) else
  setkeyv(dt.nuc, c(s.met.site, s.met.tracklabel.uni, s.met.time))

#dt.nuc.sel = myTrajExtr(in.dt = dt.nuc, in.aggr = TRUE, in.aggr.cols = c('Location_Center_X', 'Location_Center_Y'), in.met.tracklabel = s.met.tracklabel)
#dt.nuc.sel = myTrajExtr2(in.dt = dt.nuc, in.met.tracklabel = s.met.tracklabel, in.met.series = s.met.site.uni)
dt.nuc.sel = myTrajExtr3(in.dt = dt.nuc, in.met.tracklabel = s.met.tracklabel.uni)


# add real time column
#dt.nuc.sel = merge(dt.nuc.sel, dt.tass, by = s.met.time)

# add experimental description to Nuclei
if (b.well) {
  dt.nuc.sel = merge(dt.nuc.sel, dt.exp, by.x = s.met.well, by.y = 'Well')
  
} else {
  if (l.par$exp.wound)
    dt.nuc.sel = merge(dt.nuc.sel, dt.exp[, .(Position, WoundDirection, Stim_All_C, Stim_All_S)], by.x = s.met.site, by.y = 'Position') else
      dt.nuc.sel = merge(dt.nuc.sel, dt.exp[, .(Position, Stim_All_C, Stim_All_S)], by.x = s.met.site, by.y = 'Position')
}


# number of cells per site
# later saved as pdf tables
if (b.well)
dt.ncells.persite = dt.nuc.sel[, .SD[1], by = TrackObjects_Label_20_uni][, .N, by = Metadata_Well] else
  dt.ncells.persite = dt.nuc.sel[, .SD[1], by = TrackObjects_Label_20_uni][, .N, by = Position]
    
dt.ncells.percond = dt.nuc.sel[, .SD[1], by = TrackObjects_Label_20_uni][, .N, by = Condition]



####
## Calculation
# before calculating dXY/dt, it's important to order rows
if (b.well)
  setkeyv(dt.nuc.sel, c(s.met.well, s.met.site, s.met.tracklabel, s.met.time)) else
  setkeyv(dt.nuc.sel, c(s.met.site, s.met.tracklabel, s.met.time))

# calculate dXY (shift between time points)
# also, dXY is converted from pixel/l.par$exp.acquisition.freq(min) to pixel_size(um)/l.par$exp.acquisition.freq
# velocity is given in 1/h
c.scale.vel = l.par$pixel.size / l.par$exp.acquisition.freq * 60
c.scale.dist = l.par$pixel.size 

dt.nuc.sel[, dX := (-Location_Center_X + shift(Location_Center_X, 1, type = 'lead')), by = s.met.tracklabel.uni]
dt.nuc.sel[, dY := (-Location_Center_Y + shift(Location_Center_Y, 1, type = 'lead')), by = s.met.tracklabel.uni]
dt.nuc.sel[, dR := (dX^2 + dY^2)^.5]

# calcluate instantenous direction on [0, 2pi]
# Coordinate system assumes the origin [0,0] in the upper left corner
# This reflects the coordinate system of the image!
dt.nuc.sel[, dTh := myAtan(dX, dY) - pi/2, by = 1:nrow(dt.nuc.sel)]
dt.nuc.sel[, dV.umh := dR * c.scale.vel]

# median intantaneous velocity per condition
dt.nuc.sel.stat.instvel.med.percond = dt.nuc.sel[, median(dV.umh, na.rm = TRUE), by = Condition]

dt.pers = dt.nuc.sel[, .( dRdirect = ((last(Location_Center_X) - first(Location_Center_X))^2 + (last(Location_Center_Y) - first(Location_Center_Y))^2)^.5,
                dRtotal  = sum(dR, na.rm = TRUE),
                dV.umh   = mean(dV.umh, na.rm = TRUE)), 
           by = c('Condition', s.met.tracklabel.uni)][, pers := dRdirect / dRtotal]

dt.pers[, dRtotal.um := dRtotal * c.scale.dist]

  





# Shift X & Y with respect to first time point
# This is for XY plot of trajectories where all start at [0, 0]
# As above, the y-axis is mirrored w.r.t. horizontal line to reflect
# the coordinate system with [0,0] in the upper left corner.
dt.nuc.sel[, `:=`(Location_Center_X.shift = Location_Center_X - first(Location_Center_X),
                  Location_Center_Y.shift = - Location_Center_Y + first(Location_Center_Y)), by = s.met.tracklabel.uni]




if (l.par$exp.wound) {
  # rotate trajectories to follow the wound 
  # (based WoundDirection column in experimental description)
  # The wound is in S (pi) direction
  dt.nuc.sel[, c('Location_Center_X.shift', 
                 'Location_Center_Y.shift') := myRotate(Location_Center_X.shift, 
                                                            Location_Center_Y.shift, 
                                                            (- WoundDirection + 180)* pi / 180)]
}





# calculate centre of mass, COM (sum of all trajectory endpoints)
dt.com.end = rbind(
  dt.nuc.sel[, .(last(Location_Center_X.shift), 
                 last(Location_Center_Y.shift)), 
             by = c('Condition', s.met.tracklabel.uni)][, .(Location_Center_X.shift = mean(V1), 
                                                             Location_Center_Y.shift = mean(V2), 
                                                             COM = "end"), by = Condition],
  data.table(Condition = unique(dt.nuc.sel$Condition),
             Location_Center_X.shift = c(0,0),
             Location_Center_Y.shift = c(0,0),
             COM = rep('begin',2))
)




#####
## Plotting 


p.out.cond = list()
p.out.site = list()
p.out.box = list()

p.out.cond$plot_xy_perCond = myGgplotXYtraj(
  dt.nuc.sel,
  x.arg = 'Location_Center_X.shift * c.scale.dist',
  y.arg = 'Location_Center_Y.shift * c.scale.dist',
  group.arg = s.met.tracklabel.uni,
  facet.arg = 'Condition',
  facet.ncol.arg = l.par$plot.cond.facets.ncol,
  xlab.arg = '\nPosition X (um)',
  ylab.arg = 'Position Y (um)\n',
  plotlab.arg = s.exp.name, 
  path.size.arg = 0.5, 
  point.arg = dt.com.end, 
  point.col.arg = 'COM',
  point.size.arg = 3
)

p.out.cond$plot_xy_colDir_perCond = myGgplotXYtraj(
  dt.nuc.sel,
  x.arg = 'Location_Center_X * c.scale.dist',
  y.arg = 'Location_Center_Y * c.scale.dist',
  y.rev.arg = TRUE,
  group.arg = s.met.tracklabel.uni,
  facet.arg = 'Condition',
  facet.ncol.arg = l.par$plot.cond.facets.ncol,
  xlab.arg = '\nPosition X (um)',
  ylab.arg = 'Position Y (um)\n',
  plotlab.arg = s.exp.name,
  path.col.arg = "dTh",
  path.size.arg = 2
)

p.out.site$plot_xy_colDir_perSite = myGgplotXYtraj(
  dt.nuc.sel,
  x.arg = 'Location_Center_X * c.scale.dist',
  y.arg = 'Location_Center_Y * c.scale.dist',
  y.rev.arg = TRUE,
  group.arg = s.met.tracklabel.uni,
  facet.arg = s.met.well,
  facet.ncol.arg = l.par$plot.site.facets.ncol,
  xlab.arg = '\nPosition X (um)',
  ylab.arg = 'Position Y (um)\n',
  plotlab.arg = s.exp.name,
  path.col.arg = "dTh",
  path.size.arg = 2
)

p.out.cond$plot_instV_perCond = myGgplotTraj(
  dt.arg = dt.nuc.sel,
  x.arg = 'RealTime.h',
  y.arg = "dV.umh",
  group.arg = s.met.tracklabel.uni,
  facet.arg = 'Condition',
  facet.ncol.arg = l.par$plot.cond.facets.ncol,
  summary.arg = TRUE,
  xlab.arg = "\nTime (h)",
  ylab.arg = "Instantaneous velocity (um/h)\n",
  plotlab.arg = s.exp.name,
  maxrt.arg = l.par$exp.analysis,
  xaxisbreaks.arg = 2
)


p.out.box$plot_box_pers_perCond = myGgplotBoxPlot(dt.pers, 
                x.arg = 'Condition',  
                y.arg = 'pers',
                group.arg = 'Condition',
                xlab.arg = '',
                ylab.arg = 'Persistence\n',
                plotlab.arg = s.exp.name,
                y.lim.perc = c(0.05, 0.95))


p.out.box$plot_box_totDist_perCond = myGgplotBoxPlot(dt.pers, 
                x.arg = 'Condition',  
                y.arg = 'dRtotal.um',
                group.arg = 'Condition',
                xlab.arg = '',
                ylab.arg = 'Total distance (um)\n',
                plotlab.arg = s.exp.name,
                y.lim.perc = c(0.05, 0.95))


p.out.box$plot_box_instVel_perCond = myGgplotBoxPlot(dt.pers, 
                                                     x.arg = 'Condition',  
                                                     y.arg = 'dV.umh',
                                                     group.arg = 'Condition',
                                                     xlab.arg = '',
                                                     ylab.arg = 'Mean instantaneous velocity (um/h)\n',
                                                     plotlab.arg = s.exp.name,
                                                     y.lim.perc = c(0.05, 0.95))

p.out.box$plot_dens_instV_perCond = myGgplotDens(dt.nuc.sel[complete.cases(dt.nuc.sel)], 
             x.arg = 'dV.umh',  
             group.arg = 'Condition',
             xlab.arg = '\nInstantaneous velocity (um/h)',
             ylab.arg = 'Density\n',
             plotlab.arg = s.exp.name,
             vline.dt.arg = dt.nuc.sel.stat.instvel.med.percond)


# Interactive plot
# if (interactive()) {
#    library(plotly)
# 
#   (gg <- ggplotly(p.out.cond$plot_xy_perCond))
#   (gg <- ggplotly(p.out.cond$plot_xy_colDir_perCond))
#   (gg <- ggplotly(p.out.site$plot_xy_colDir_perSite))
#   (gg <- ggplotly(p.out.cond$plot_instV_perCond))
# }
 # 
 # 
 # 
#####
## Save files

# Create directory for plots in the currenty working directory
ifelse(!dir.exists(file.path(".", l.par$dir.plot)), dir.create(file.path(".", l.par$dir.plot)), FALSE)

if (l.par$plot.save) {
  # all plots
  lapply(names(p.out.cond),
         function(x)
           ggsave(
             filename = paste0(l.par$dir.plot, '/', x, ".pdf"),
             plot = p.out.cond[[x]],
             width = l.par$plot.cond.width, 
             height = l.par$plot.cond.height
           ))
  
  lapply(names(p.out.site),
         function(x)
           ggsave(
             filename = paste0(l.par$dir.plot, '/', x, ".pdf"),
             plot = p.out.site[[x]],
             width = l.par$plot.site.width, 
             height = l.par$plot.site.height
           ))
  
  
  lapply(names(p.out.box),
         function(x)
           ggsave(
             filename = paste0(l.par$dir.plot, '/', x, ".pdf"),
             plot = p.out.box[[x]],
             width =  l.par$plot.box.width, 
             height = l.par$plot.box.height
           ))

    # tables
  # tab-delimited files per condition that are read by IBIDI tool
  dt.nuc.out = dt.nuc.sel[, c('Condition', s.met.tracklabel.uni, 'RealTime.m', 'Location_Center_X', 'Location_Center_Y'), with = FALSE]
  dt.nuc.out[, `:=`(rowIndx = .I, TrackNo = .GRP, SliceNo = 1:.N, PosX = as.integer(round(1000*Location_Center_X)), PosY = as.integer(round(1000*Location_Center_Y))), by = s.met.tracklabel.uni]
  
  lapply(unique(dt.nuc.out[, Condition]),
         function(x)
           write.table(file = paste0(l.par$dir.plot, '/', "tab_tracks_", gsub(" ", "", x), ".txt"), x = dt.nuc.out[Condition == x, .(rowIndx, TrackNo, SliceNo, PosX, PosY)], row.names = FALSE, quote = FALSE, sep = "\t", eol = "\r\n")
         )
    
  # number of cells per site and per condition
  pdf(paste0(l.par$dir.plot, '/', "tab_ncells_perSite.pdf"),
      height = 7,
      width = 5)
  grid.table(dt.ncells.persite)
  dev.off()
  
  pdf(paste0(l.par$dir.plot, '/', "tab_ncells_perCondition.pdf"),
      height = 4,
      width = 4)
  grid.table(dt.ncells.percond)
  dev.off()
}

