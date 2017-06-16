# checks subpopulations based on bimodal distribution
# of instantaneous velocities

# execute main script first!

myGgplotDens(dt.nuc.sel[complete.cases(dt.nuc.sel)], 
             x.arg = 'dV.umh',  
             group.arg = 'Condition',
             xlab.arg = '\nInstantaneous velocity (um/h)',
             ylab.arg = 'Density\n',
             plotlab.arg = s.exp.name,
             vline.dt.arg = dt.nuc.sel.stat.instvel.med.percond)


nrow(dt.nuc.sel)
nrow(dt.nuc.sel[complete.cases(dt.nuc.sel)])

dt.nuc.sel.comp = dt.nuc.sel[complete.cases(dt.nuc.sel)]

dt.nuc.sel.sub = dt.nuc.sel.comp[Condition %like% 'Chemokinesis']

myGgplotDens(dt.nuc.sel.sub, 
             x.arg = 'dV.umh',  
             group.arg = 'Condition',
             xlab.arg = '\nInstantaneous velocity (um/h)',
             ylab.arg = 'Density\n',
             plotlab.arg = s.exp.name,
             vline.dt.arg = dt.nuc.sel.stat.instvel.med.percond)

dt.nuc.sel.sub[, peak := ifelse(dV.umh < 7, 0, 1)]

dt.nuc.sel.sub[, mean(peak), by = 'Condition']
