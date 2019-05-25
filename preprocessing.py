#!/usr/bin/env python
import numpy
import math
import random
#import matplotlib.pyplot as plt
import fileinput
import os
import json
import sys
import re

class Atom:
    def __init__(self, line, move, theta):
        self.serial = int(line[6:11])
        self.name = line[11:16].strip()
        self.resName = line[17:20]
        self.resID = int(line[22:26])
        self.x = float(line[27:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])
        self.occupancy = line[55:60].strip()
        #self.radius = float(line[60:66].strip())
        if self.occupancy: self.occupancy = float(self.occupancy)

        #move and rotate
        rotate = self.point_transfer(self.x, self.y, self.z, move, theta)
        self.x = rotate[0, 0]
        self.y = rotate[1, 0]
        self.z = rotate[2, 0]

    def point_transfer(self, x, y, z, move, theta):
        new_point = numpy.zeros((3, 1), dtype = float)
        #move
        new_point[0, 0] = x + move
        new_point[1, 0] = y + move
        new_point[2, 0] = z + move
        #rotate
        R_x = numpy.array([[1, 0, 0],[0, math.cos(theta[0]), -math.sin(theta[0])],
                    [0, math.sin(theta[0]), math.cos(theta[0])]])           
        R_y = numpy.array([[math.cos(theta[1]), 0, math.sin(theta[1])], [0, 1, 0 ],
                    [-math.sin(theta[1]), 0, math.cos(theta[1])]])
        R_z = numpy.array([[math.cos(theta[2]), -math.sin(theta[2]), 0],
                    [math.sin(theta[2]), math.cos(theta[2]), 0], [0, 0, 1]])
        res = numpy.dot(R_z, numpy.dot(R_y, numpy.dot(R_x, new_point)))
        return res

    def __getitem__(self, key):
        return self.__dict__[key]

class PDB:
    def __init__(self, file, h):
        self.h = h
        self.file = file
        self.atoms = []
        self.central = []
        self.range = []
        self.gridpointnum = [120, 120, 120]
        self.terminal = 0
        self.parse()
        self.get_range()
        self.get_central()
        

    def parse(self):
        MODEL = None
        move = random.uniform(-10.0, 10.0)
        theta = [random.randint(0, 360), random.randint(0, 360), random.randint(0, 360)]
        print ("move: ", move, "theta: ", theta)
        f = open(self.file, 'r')
        for line in f.readlines():
            if line.startswith('ATOM'):
                atom = Atom(line, move, theta)
                self.atoms.append(atom)
        f.close()

        f = open(self.file, 'r')
        lines = f.readlines()
        print ("filename", self.file)
        try:
            self.terminal = int(lines[len(lines) - 4][23:26])
        except ValueError:
            self.terminal = int(lines[len(lines) - 10][23:26])
            
        f.close()
        
        self.overwrite()

    def getterminal(self):
        return self.terminal

    def overwrite(self):
        i = 0
        for line in fileinput.FileInput(self.file, inplace = 1): 
            if line.startswith('ATOM'):
                tmp_x = ("%.3f" % self.atoms[i].x)
                tmp_y = ("%.3f" % self.atoms[i].y)
                tmp_z = ("%.3f" % self.atoms[i].z)
                print (line.replace(line[29:38],str(tmp_x).rjust(9," ")).replace(line[38:46],str(tmp_y).rjust(8," ")).replace(line[46:54],str(" "+tmp_z).rjust(8," ")).rstrip('\n'))
                i += 1
            else:
                print (line.rstrip('\n'))

    def get_atoms(self, to_dict=True):
        """Return a list of all atoms.

        If to_dict is True, each atom is represented as a dictionary.
        Otherwise, a list of Atom objects is returned."""
        if to_dict: return [x.__dict__ for x in self.atoms]
        else: return self.atoms

    def get_model(self, model_num, to_dict=True):
        """Return all atoms where MODEL == model_num"""
        model_atoms = [x for x in self.atoms if x.MODEL == model_num]
        if to_dict:
            return [atom.__dict__ for atom in model_atoms]
        else:
            return model_atoms
    
    def get_range(self):
        # get maximum/minimum coordinate value of X,Y,Z
        x = []
        y = []
        z = []
        for atom in self.atoms:
            x.append(atom.x)
            y.append(atom.y)
            z.append(atom.z)
        self.range =  [max(x), min(x), max(y), min(y), max(z), min(z)]
             
    def get_central(self):
        #calculate center of grid
        ran = self.range
        self.central = [(ran[0] + ran[1])/2, (ran[2] + ran[3])/2, (ran[4] + ran[5])/2]
        
    
    def get_index(self, position):
        # get coordinate position on grid  
        N = self.gridpointnum
        central = self.central
        index = []
        for i in range(0,3):
            index.append(int(math.floor((position[i] - central[i] + 60*self.h)/self.h)))
        return index
    
    def get_position(self, index):
        #get coordinate position by given index
        N = self.gridpointnum
        central = self.central
        position = []
        for i in range(0,3):
            position.append(central[i] - 60*self.h + index[i]*self.h)
        return position


    def get_grid(self, atoms):
        # To determine whether a point included in a atom
        i = int(self.gridpointnum[0])
        j = int(self.gridpointnum[1])
        k = int(self.gridpointnum[2])
        grid = numpy.zeros((i,j,k), dtype = numpy.float64)
        cal = 0
        vol = 0
        #print 'checking..'
        for atom in atoms:
            #vol = vol + 4.18*atom.radius*atom.radius*atom.radius
            #print("check1")
            #build the box
            lower_position = [atom.x - 2, atom.y - 2, atom.z - 2]
            higher_position = [atom.x + 2, atom.y + 2, atom.z + 2]
            #print("check2")
            lower_index = self.get_index(lower_position)
            higher_index = self.get_index(higher_position)
            #print("check3")
            flag = False
            for pos in lower_index:
                if (pos >= 60) or (pos <= -60):
                    flag = True
            for pos in higher_index:
                if (pos >= 60) or (pos <= -60):
                    flag = True
            if flag == True:
                continue
            #print("check4")
            for x in range(lower_index[0], higher_index[0]):
                for y in range(lower_index[1], higher_index[1]):
                    for z in range(lower_index[2], higher_index[2]):
                        if(grid[x, y, z] != 0):
                                continue
                        #print("check5")
                        #distance
                        position = self.get_position([x, y, z])
                        #print position, atom.x, atom.y, atom.z, "tmp"
                        dist = numpy.sqrt(pow((atom.x - position[0]),2) + pow((atom.y - position[1]),2) + pow((atom.z - position[2]),2))
                        #print dist, atom.radius, "radius"
                        #print("check6")
                        if dist <= 2.0:
                            pr = math.exp(-(dist*dist/2))
                            grid[x,y,z] += pr

        #print 'checking complete'
        #print vol, "volume"
        #print("check7")
        return grid

    def volume_cal(self, atoms):
        result = 0

        i = int(self.gridpointnum[0])
        j = int(self.gridpointnum[1])
        k = int(self.gridpointnum[2])
        grid = numpy.zeros((i,j,k), dtype = numpy.float64)

        grid = self.get_grid(atoms)

        for x in range(0, i):
            for y in range(0, j):
                for z in range(0, k):
                    result = result + grid[x, y, z]
        result = result * pow(self.h, 3)
        return result


def classify_atom(atom, terminal_resid):
    if (atom.resName == 'CYG' and atom.name == 'SG') or (atom.resName == 'MET' and atom.name == 'SD') \
    or (atom.resName == 'MSE' and atom.name == 'SE'):
        return 0
    if (atom.resName == 'ASN' and atom.name == 'ND2') or (atom.resName == 'GLN' and atom.name == 'NE2') or (atom.name == 'N'):
        return 1
    if (atom.resName == 'HIS' and (atom.name == 'ND1' or atom.name == 'NE1')) or (atom.resName == 'TRP' and atom.name == 'NE1'):
        return 2
    if atom.resName == 'ARG' and (atom.name == 'NE' or atom.name.startswith('NH') == True):
        return 3
    if atom.resName == 'LYS' and atom.name == 'NZ':
        return 4
    if (atom.resName == 'ASN' and atom.name == 'OD1') or (atom.resName == 'GLN' and atom.name == 'OE1') or (atom.name == 'O' and atom.resID != terminal_resid):
        return 5
    if (atom.resName == 'SER' and atom.name == 'OG') or (atom.resName == 'THR' and atom.name == 'OG1') or \
    (atom.resName == 'TYR' and atom.name == 'OH'):
        return 6
    if (atom.resName == 'ASP' and atom.name.startswith('OD') == True) or (atom.resName == 'GLU' and atom.name.startswith('OE') == True) or\
    (atom.resID == terminal_resid and atom.name == 'O') or (atom.resID == terminal_resid and atom.name == 'OXT'):
        return 7
    if (atom.resName == 'ARG' and atom.name == 'CZ') or (atom.resName == 'ASN' and atom.name == 'CG') \
    or (atom.resName == 'ASP' and atom.name == 'CG') or (atom.resName == 'GLN' and atom.name == 'CD') or (atom.resName == 'GLU' and atom.name == 'CD') or (atom.name == 'C'):
        return 8
    if (atom.resName == 'HIS' and (atom.name == 'CG' or atom.name == 'CD2' or atom.name == 'CE1'))  \
    or (atom.resName == 'PHE' and (atom.name.startswith('CD') == True or atom.name.startswith('CE') == True or atom.name == 'CZ' or atom.name == 'CG')) \
    or (atom.resName == 'TRP' and (atom.name.startswith('CD') == True or atom.name.startswith('CE') == True or atom.name.startswith('CE') == True  or atom.name == 'CG' or atom.name == 'CH2'))\
    or (atom.resName == 'TYR' and (atom.name.startswith('CD') == True or atom.name.startswith('CE') == True or atom.name == 'CZ' or atom.name == 'CG')) :
        return 9
    if (atom.resName == 'ALA' and atom.name == 'CB') or (atom.resName == 'ARG' and (atom.name == 'CB' or atom.name == 'CG' or atom.name == 'CD'))\
    or ((atom.resName == 'ASN' or atom.resName == 'ASP' or atom.resName == 'CYS') and atom.name == 'CB')\
    or (atom.resName == 'GLN' and (atom.name == 'CB' or atom.name == 'CG')) or (atom.resName == ' GLU' and (atom.name == 'CB' or atom.name == 'CG')) \
    or (atom.resName == 'HIS' and atom.name == 'CB') or (atom.resName == 'ILE' and (atom.name == 'CB' or atom.name.startswith('CG') or atom.name == 'CD1'))\
    or (atom.resName == 'LEU' and (atom.name == 'CB' or atom.name == 'CG' or atom.name.startswith('CD'))) \
    or (atom.resName == 'LYS' and (atom.name == 'CB' or atom.name == 'CG' or atom.name == 'CD' or atom.name == 'CE')) \
    or (atom.resName == 'MET' and (atom.name == 'CB' or atom.name == 'CG' or atom.name == 'CE'))\
    or (atom.resName == 'MSE' and (atom.name == 'CB' or atom.name == 'CG' or atom.name == 'CE'))\
    or ((atom.resName == 'PHE' or atom.resName == 'SER' or atom.resName == 'TRP' or atom.resName == 'TYR') and atom.name == 'CB') \
    or (atom.resName == 'PRO' and (atom.name == 'CB' or atom.name == 'CG' or atom.name == 'CD'))\
    or (atom.resName == 'THR' and (atom.name == 'CB' or atom.name == 'CG2')) or (atom.resName == 'VAL' and (atom.name == 'CB' or atom.name.startswith('CG'))) or (atom.name == 'CA'):
        return 10
    return -1

def getDecoysList(decoysPath):
    decoysInfo = decoysPath + "/list.dat"
    f = open(decoysInfo, 'r')
    min_gdtts = 1000
    max_gdtts = -1000
    for line in f.readlines()[1:]:
        data = line.split("\t")
        gdtts = float(data[3].strip())
        if gdtts < min_gdtts:
            min_gdtts = gdtts
        if gdtts > max_gdtts:
            max_gdtts = gdtts
    f.close()
    
    bins = {}
    res = []
    Nb = 10

    f = open(decoysInfo, 'r')
    for line in f.readlines()[1:]:
        data = line.split("\t")
        gdtts = float(data[3].strip())
        n = 1 + math.floor(Nb * (gdtts - min_gdtts) / (max_gdtts - min_gdtts))
        if n not in bins:
            bins[n] = []
            bins[n].append(data[0].strip())
        else:
            bins[n].append(data[0].strip())

    for key in bins:
        random.shuffle(bins[key])
        tmp = decoysPath + "/" + bins[key][0]
        res.append(tmp)
    print ("number of selected decoys in this target : ", len(res))
    return res

def get_featureandlable(file):
    res = []
    
    try:
        matchObj = re.match(r'(.*)/T([0-9]{4}/)(.*)', file)
        if matchObj:
            filename = matchObj.group(3)
            rootname = file[:46] + matchObj.group(1) +"/T"+ matchObj.group(2) + "list.dat"
        
        gdtts = -1
        f = open(rootname, 'r')
        for line in f.readlines()[1:]:
            data = line.split("\t")
            if data[0] == filename:
                gdtts = float(data[3].strip())
                break
                
        res.append(gdtts)        
        #return res
        return gdtts
        
    except:
        print ("pass a value exception")
        res.append(0.0)
        #return res
        return 0.0


def getFilename(targetpath, target_num):
    all_filenames = []
    decoy_nums = []
    
    targetPath = targetpath 
    targetList=list(os.path.join(targetPath,name) for name in os.listdir(targetPath))
    randomTarget=list(random.sample(targetList, target_num)) #change this number to get more targets.

    for decoysPath in randomTarget:
        randomDecoy = getDecoysList(decoysPath)
        all_filenames = all_filenames + randomDecoy
        decoy_nums.append(len(randomDecoy))
        
    #train_filenames, test_filenames = split_trainingANDtesting(all_filenames)

    #return train_filenames, test_filenames
    return all_filenames, decoy_nums


def getData(filename):
    
    file = filename.decode('ascii')
    file = str(file)
    
    #pdb = PDB(file, float(sys.argv[1]))
    pdb = PDB(file, 1.0) # h = 1
    atoms = pdb.get_atoms(to_dict=False)
    terminal_resid = pdb.getterminal()

    atom_all = []
    print ("range: ", pdb.range)
    print ("central: ", pdb.central)
    print ("gird number: ", pdb.gridpointnum[0], pdb.gridpointnum[1], pdb.gridpointnum[2])
     # print the name of each atom

    classified_atoms = [[],[],[],[],[],[],[],[],[],[],[]]
    for atom in atoms:
        atom_all.append(atom)
        n = classify_atom(atom, terminal_resid)
        if n == -1:
            continue
        else:
            classified_atoms[n].append(atom)

    feature = []
    for i, atom_list in enumerate(classified_atoms):
        tmp = pdb.get_grid(atom_list)
        feature.append(tmp.tolist())

    print ("number of atoms: ", len(atom_all))
    print (numpy.shape(feature))
    #print 'volume : ', pdb.volume_cal(atom_all)
    gdtts = get_featureandlable(file)
    
    one_hot = numpy.zeros(10,dtype = numpy.int32)
    index = int(math.ceil(gdtts/0.1))
    one_hot[index - 1] = 1
    
    return feature, one_hot