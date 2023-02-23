## THIS CODE RUNS THE AGENT-BASED, STOCHASTIC MODEL FOR A STARTING NUMBER OF ENGRAFTED CELLS (CALCULATED IN THE CODE AS DOSE X CLONE.NO)
## THAT GO ON TO PROLIFERATE AND DIFFERENTIATE FOR EACH RANDOM SEED ("SEED.NO") OVER THE TIME SPAN SPECIFIED IN 
## "RUN.TIME." "CYCLE.PERC" IS THE PERCENTAGE OF CELLS THAT DIVIDE AT EACH TIME STEP. THE MODEL TRACKS THE 
## TOTAL NUMBER OF CELLS IN THE VECTOR "WORK.CELLS" THE MODEL ALLOWS EACH CELL THAT DIVIDES INTO TWO COPIES TO 
## AT EACH TIME STEP TO EITHER ADD AN IDENTICAL CELL (PROLIFERATION) OR SUBTRACTS A CELL (DIFFERENTIATION) WITH 
## 50% PROBABILITY ON AVAERAGE. BECAUSE BOTH COPIES OF A DIVIDED CELL CAN HAVE EITHER FATE, THIS ALLOWS FOR 
## SYMMETRIC DIVISION.

## THE DATA FRAME "FREQ" STORES THE TOTAL NUMBER OF CELLS AND THE TOTAL NUMBER OF CLONES OF EACH TYPE 
## FOR EACH RANDOM SEED

## IN THIS SCRIPT, THE FRACTION (CYCLE.PERC) OF CELLS CYCLING AT EACH TIME STEP IS 0.15; THERE ARE 4 STARTING ENGRAFTED CELL
## NUMBERS (100K, 200K, 400K, 600K). 100 SIMULATIONS ARE RUN AND PLOTTED FOR EACH INITIAL CONDITION.

library(scales)
library(dplyr)
library('RColorBrewer')
library('lubridate')
library('ggplot2')
library(wesanderson)
library(ggpubr)

rm(list=ls())
graphics.off()

clone.no = 10^5
run.time = 1000
seed.no = 100
cycle.perc = 0.15
doses = 4

txplt = c(1, 2, 4, 6) 

fig.freq.png = paste0("20230208_multi_dose_seeds_",  clone.no, "_clones_",  run.time, "_days_", seed.no,  "_seeds_", cycle.perc, "_cycling", ".png")
fig.freq.svg = paste0("20230208_LINEAR_multi_dose_",  clone.no, "_clones_",   seed.no, "_seed_", run.time, "_days", cycle.perc, "_cycling",".svg")
csv.name = paste0("20230208_multi_dose_", clone.no, "_clones_", run.time, "_days_", seed.no,  "_seeds_", cycle.perc, "_cycling",".csv")

#clone storage data frame
days = rep(1:run.time, seed.no*doses)
clones = rep(0, seed.no*run.time*doses)
log.clones = rep(0, seed.no*run.time*doses)
cells = rep(0, seed.no*run.time*doses)
log.cells = rep(0, seed.no*run.time*doses)
seed = rep(0, seed.no*run.time*doses)
dose = rep(0, seed.no*run.time*doses)

freq = data.frame(days, dose, seed, clones, cells, log.clones, log.cells)

doser = c()  #determine cell doses to pre-fill in dataframe
for (m in 1:doses){
  doser = c(doser, rep(txplt[m]*clone.no, run.time*seed.no))
}

seeder = c() #determine seeds to pre-fill in df
for (h in 1:seed.no){
  seeder = c(seeder, rep(h, run.time))
}

#add to df
freq$dose = doser 
freq$seed = rep(seeder, doses)

start = Sys.time() #time to run 

for (i in 1:doses){ #for each dose
  dosee = txplt[i]
  
  #run the model for each seed.no random seeds
  for (j in 1:seed.no){ 
    seeder = j
    set.seed(seeder)
  
    work.cells = 1:(dosee*clone.no) #total starting cell number
    diffed.cells = c()
  
    for (k in 1:(run.time)){ #for each time in run.time
      # choose cells to divide
      div.index = sample(1:length(work.cells), floor(length(work.cells) * cycle.perc),replace = FALSE)
      div.cells = work.cells[div.index]
      work.cells = work.cells[-div.index] #take the cycling cells out of the work.cells vector
      # duplicate div cells
      doubled = c(div.cells,div.cells)
      # determine fate - stay stem or differentiate
      fate = runif(length(doubled))
      new.stem.index = which(fate >= 0.5)
      new.diff.index = which(fate < 0.5)
      # put proliferated and differentiated cells in vectors
      new.stem = doubled[new.stem.index]
      new.diff = doubled[new.diff.index]
      # save differentiated cells in vector with prior differentiated cells
      diffed.cells = c(diffed.cells, unique(new.diff))
      # put the cycling cells that renewed back in the vector
      work.cells = c(work.cells, new.stem)
    
      unique.cells = length(unique(work.cells)) # count remaining unique clones, N_t
      total.cells = length(work.cells) # count remaining total cells, H_t
    
    freq$clones[((j-1)*run.time + run.time*seed.no*(i-1) + k)] = unique.cells #write these into the data frame on the next row
    freq$cells[((j-1)*run.time + run.time*seed.no*(i-1) + k)] = total.cells 
    
    }
  }
}

freq$log.clones= ifelse(freq$clones == 0, 0, log10(freq$clones))
freq$log.cells = ifelse(freq$cells == 0, 0, log10(freq$cells))

stop = Sys.time() #time to run 
stop - start



##===================LINEAR SCALE -= LINE PLOT OF UNIQUE CLONES AND TOTAL CELLS VERSUS TIME
#
freq$seed = as.factor(freq$seed)
freq$dose = as.factor(freq$dose)

ggplot(data = freq) +
  geom_line(aes(x = days, y = clones, group = interaction(dose, seed), color = dose, linetype = seed), size = 0.5) +
  coord_cartesian(ylim = c(0, 100000)) +
  xlab("Time (days)") + ylab("Number of unique clones and Total Cells") +
  scale_x_continuous(breaks = c(0, 200, 400, 600, 800, 1000)) +
  scale_y_continuous(breaks = seq(0, 100000, 20000), labels = c("0", "20000", "40000", "60000", "80000", "100000")) +
  scale_color_discrete("Starting cell number", labels = c("100K", "200K", "400K", "600K")) +
  scale_linetype_discrete(guide = "none") + 
  xlab("Days after engraftment") + ylab("Number of unique clones") +
  theme_classic() +
  theme(legend.position = c(0.6, 0.8), legend.text = element_text(size = 12), legend.title = element_text(size = 12),
        legend.key.width = unit(1, "cm"),
        axis.text = element_text(size = 12), axis.title = element_text(size = 12))

linear.plot.pdf = paste0("20230208_LINEAR_multi_dose_seeds_",  clone.no, "_clones_",  run.time, "_days_", seed.no,  "_seeds_", cycle.perc, "_cycling", ".pdf")

ggsave(linear.plot.pdf, width = 7, height = 5, units = "in", dpi = 300)

