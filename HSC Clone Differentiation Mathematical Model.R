#
#  very simple simulation for stem cells in NHP translpant
#
#  from RNAseq data:
#   snap-shot stem cell cell division rate ~= 0.15
#   all stem cells are equivalent
#   many stem cell (ISA cell ID) are lost early
#   stem cells (ISA cell ID) become stabily identified late
#    
#   to simulate will assume 1 day = amount of time for stem cell to complete cell cycle
#   all cells equally likely to cycle
#   when a cell cycles, each daughter makes a decision to stay stem or differentiate
#   odds are 50-50  for a stable stem cell population
#
#
# running this script shows stem cell diversity drops from 10,000 to about 700 
# with some small run-to-run variation 
#
# the script shows large clones of about 100 members and lower form
#
# the script shows many clones fully differentiate early and are then lost
#   
#
#
INIT_STEM = 120000

{
  cells = 1:INIT_STEM
  all.dif.cells = list()
  clones = c()
  
  for (day in 1:365) {
    # pick cells to divide
    div.index = sample(1:length(cells),floor(length(cells) * 0.15),replace = FALSE)
    div.cells = cells[div.index]
    cells = cells[-div.index]
    # duplicate div cells
    ly = c(div.cells,div.cells)
    # fate, stay stem or differentiate
    fate = runif(length(ly))
    new.stem.index = which(fate >= 0.5)
    new.diff.index = which(fate < 0.5)
    new.stem = ly[new.stem.index]
    new.diff = ly[new.diff.index]
    # save diff cells (collapse multiple)
    all.dif.cells = append(all.dif.cells,list(unique(new.diff)))
    # add new stem cells
    cells = c(cells,new.stem)
    
    unique.cells = length(table(cells))
    clones = c(clones,unique.cells)
    
    print(paste0("day ",day," number of stem cells = ",length(cells)," u = ",unique.cells," dif = ",length(new.diff)))
    
  }
}
#
# unique clones by day
#
#
plot(clones,type='l',ylab = "Unique Clones",xlab = "Day",main="Unique Clones")
#
# clone sizes
#
length(cells)
length(table(cells))
r = rle(sort(cells))
plot(sort(r$lengths,decreasing = TRUE),type='l',ylab = "Clone Sizes",xlab = "Clone Index",main = "Clone Size Distribution")
sort(r$lengths,decreasing = TRUE)[1:10]
#
# for each day, how many cells are seen one time
#
final = unlist(all.dif.cells)
# find out how many times each cell forms a differentiated cell
# id = number of times a stem made 1 or more diff cells in a day
id = rep(0,INIT_STEM)
for (cell.id in final) {
  id[cell.id] = id[cell.id] + 1
}
sort(id)[1:10]
#
# how many times was each diff cell seen
#
lt = lapply(1:365,function(x) {
  l = unique(all.dif.cells[[x]])
  return(length(which(id[l] <= 10))/length(l))
})
plot(unlist(lt),type='l',ylab = "Fraction of Clones Seen One Time",xlab = "Day",main="Fraction Of Single Time Clones")
