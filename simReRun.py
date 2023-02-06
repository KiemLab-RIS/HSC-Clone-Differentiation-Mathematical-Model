#
# simulate clone development and re-run multiple times with different seeds
#
#
#  Mark Enstrom
#  1-25-23
#
#  Assume only LT HSC
#  Assume 50% HSC/Diff rate
#  Assume 15% cell cycle per day
#
#
import os
import sys
import random
from collections import defaultdict
from collections import namedtuple
from enum import Enum
#
# useage simReRun.py rep% stem%
#
if (len(sys.argv) != 3):
  print("useage simDie.py rep%  stem%")
  quit()
LT_REP  = float(sys.argv[1])
LT_STEM = float(sys.argv[2])

print(f"simReRun.py rep rate = {LT_REP}%  stem rate = {LT_STEM}%")
#
#
# stem cells have 1 type....HSC
#
# assumptions:
#    what % LT cycle
#    what % LT daughter -> STEM vs PROGENITOR
#    what % LT die after division
#    how long to sim
#    starting cells (one cell per clone)
#
#LT_REP         = 15  # percent
#LT_STEM        = 50  # percent
SIM_DAYS       = 365
STARTING_CELLS = 120000
#
# track cells
#
globalCellID = 0
#
# diversity over time and progenotors created
#
diversity = []
progenitors = []
totalStem = []
#
# cell types
#
class CellType(Enum):
    LT = 1
    PG = 3
    E  = 4
    M  = 5
    L  = 6
#-----------------------------------------
#
# cell class
#   cell has a clone id
#   cell has a type
#   cell has a replication rate....which may be based on type or individual
#       ie: faster if division count is higher
#   cell has differention rate....type or individual
#   cell has expansion...current expansion level
#   cell has expansion target
#
#-----------------------------------------
class Cell:
  classCellID = 0
  def __init__(self,cloneID,division,repRate,difRate,e1,e2,type):
    self._cloneID = cloneID
    self._cellID = Cell.classCellID
    Cell.classCellID += 1
    self._division = division
    self._repRate = repRate
    self._difRate = difRate
    self._expansion = e1
    self._expansionTarget = e2
    self._type = type
  def clone(self):
    c2 = Cell(self._cloneID,self._division,self._repRate,self._difRate,self._expansion,self._expansionTarget,self._type)
    return c2
#-----------------------------------------
#
# Bone Marrow (MB) class contains all cells
#
#-----------------------------------------
class BM:
  def __init__(self,numLT):
    # all cells
    self._cells = {}
    #
    # temp cell groups
    #
    self._simRep = {}
    self._simNoRep = {}
    self._simNewHSC = {}
    #
    #
    #
    for i in range(0,numLT):
      c1 = Cell(i,0,LT_REP,LT_STEM,0,0,CellType.LT)
      self._cells[c1._cellID] = c1
  #
  # simulate replication, all cells either end up in
  # cell-cycle list or not-cycle list
  #
  def simRep(self):
      for c1 in self._cells.values():
        r1 = random.random() * 100
        if r1 <= c1._repRate:
          c1._division += 1
          c2 = c1.clone()
          self._simRep[c1._cellID] = c1
          self._simRep[c2._cellID] = c2
        else:
          self._simNoRep[c1._cellID] = c1
      return(len(self._simRep))
  #
  # for all cell that just replicated:
  # determine if they stay the same or differentiate
  #
  def simDif(self):
      for c1 in self._simRep.values():
          r1 = random.random() * 100
          if r1 > c1._difRate:
            # cell diff to next type
            if c1._type == CellType.LT:
              c1._type = CellType.PG
            elif c1._type == CellType.PG:
              pass
            else:
              print('CellType Error')
            #
            # always add cell to new list
            #
            self._simNoRep[c1._cellID] = c1
          else:
            #
            # daughter stays HSC
            #
            self._simNewHSC[c1._cellID] = c1
  #
  # finish sim
  #    remove PG cells and record how many created
  #    move rest back to _cells list
  #
  def simFinish(self):
    numPG = 0
    self._cells = {}
    # add noRep cells back to main list
    for c1 in self._simNoRep.values():
        if c1._type != CellType.PG:
            self._cells[c1._cellID] = c1
        else:
            numPG += 1
    #
    # for rep HSC, chance they might die
    #
    for c1 in self._simNewHSC.values():
      self._cells[c1._cellID] = c1
    #
    # clear for next round
    #
    self._simNoRep = {}
    self._simRep = {}
    self._simNewHSC = {}
    #
    # calc clonal diversity
    # save diversity and progenitors created
    #
    clones = defaultdict(int)
    for c1 in self._cells.values():
      clones[c1._cloneID] += 1
    diversity.append(len(clones))
    progenitors.append(numPG)
    totalStem.append(len(self._cells))
    #
    return numPG
  #
  # simulate 1 day
  #
  def sim(self,day):
    numRep = bm.simRep()
    bm.simDif()
    pg = bm.simFinish()
    lt = self.stats()
    print(f'day {day} lt = {lt:5d} rep = {numRep} pg created = {pg}')
    return((lt,pg))
  #
  # stats
  #
  def stats(self):
    lt = 0
    for cell in bm._cells.values():
      if cell._type == CellType.LT:
          lt += 1
    return(lt)
#
# sim   10 times
#
seeds = [ 2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
31,37, 41, 43, 47, 53, 59, 61, 67, 71]

with open('simOutReRun/simReRun.tsv','w') as fr:
  for run in range(0,len(seeds)):
    random.seed(seeds[run])
    print(f'seed = {seeds[run]}')
    diversity = []
    progenitors = []
    totalStem = []
    #
    # init BM
    #
    bm = BM(STARTING_CELLS)
    fr.write('day\tlt\tpg\tseeds[run]\n')
    for day in range(1,SIM_DAYS):
      lt,pg = bm.sim(day)
      fr.write(f'{day}\t{lt}\t{pg}\n')
    #-----------------------------------------
    #
    # simulation done
    #   save diversity and progenitors
    #   calc and save clone stats
    #
    #-----------------------------------------
    cloneLT = defaultdict(int)
    for cell in bm._cells.values():
        if cell._type == CellType.LT:
            cloneLT[cell._cloneID] += 1

    sd = sorted(cloneLT.items(),key=lambda x:x[1],reverse=True)

    fileName = f"simOutReRun/simReRun_diversity_{LT_REP}_{LT_STEM}_R{run}.tsv"
    with open(fileName,'w') as fh:
        fh.write(f"# simDie.py\n")
        fh.write(f"# HSC cycle rate = {LT_REP}\n")
        fh.write(f"# HSC Stem rate = {LT_STEM}\n")
        fh.write(f"# days simulate = {SIM_DAYS}\n")
        fh.write(f"# starting cells = starting clones = {STARTING_CELLS}\n")
        fh.write(f"# POST Simulation:\n")
        fh.write(f"# total HSC    = {len(bm._cells)}\n")
        fh.write(f"# total clones = {len(sd)}\n")
        fh.write(f"Diversity\tProgenitors\n")
        for i in range(0,len(diversity)):
            fh.write(f"{diversity[i]}\t{progenitors[i]}\n")


    fileName = f"simOutReRun/simReRun_clones_{LT_REP}_{LT_STEM}_R{run}.tsv"
    with open(fileName,'w') as fh:
        fh.write(f"# simDie.py\n")
        fh.write(f"# HSC cycle rate = {LT_REP}\n")
        fh.write(f"# HSC Stem rate = {LT_STEM}\n")
        fh.write(f"# days simulate = {SIM_DAYS}\n")
        fh.write(f"# starting cells = starting clones = {STARTING_CELLS}\n")
        fh.write(f"# POST Simulation:\n")
        fh.write(f"# total HSC    = {len(bm._cells)}\n")
        fh.write(f"# total clones = {len(sd)}\n")
        fh.write(f"CloneID\tNumber\n")
        for k,v in sd:
          fh.write(f"{k}\t{v}\n")


    fileName = f"simOutReRun/simReRun_stem_{LT_REP}_{LT_STEM}_R{run}.tsv"
    with open(fileName,'w') as fh:
        fh.write(f"# simDie.py\n")
        fh.write(f"# HSC cycle rate = {LT_REP}\n")
        fh.write(f"# HSC Stem rate = {LT_STEM}\n")
        fh.write(f"# days simulate = {SIM_DAYS}\n")
        fh.write(f"# starting cells = starting clones = {STARTING_CELLS}\n")
        fh.write(f"# POST Simulation:\n")
        fh.write(f"# total HSC    = {len(bm._cells)}\n")
        fh.write(f"# total clones = {len(sd)}\n")
        fh.write(f"CloneID\tNumber\n")
        for s in totalStem:
          fh.write(f"{s}\n")

