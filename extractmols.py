import os
import pandas as pd


file = 'unalignedmols.txt'
bnxfile = 'Fico3_DLE-AcrCen_Y-sat-_7__RawMolecules.bnx'
bnxdir = 'Fico3_unalignedmols'
with open(file, 'r') as i:
    molIDS = i.readlines()
    molIDS = [int(id.strip('\n')) for id in molIDS]
for molecule in molIDS:
    cmd = "grep $'^0\t" + str(molecule) + "\t' -A 6 " + bnxfile + " > " + bnxdir + '/' + str(molecule) + '.bnx'
    # print(cmd)
    os.system(cmd)