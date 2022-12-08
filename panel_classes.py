import numpy as np
from copy import deepcopy
from collections import defaultdict

class Fluor:
    def __init__(self, name : str ,brightness: float, \
                 all_bonds = None, channel_values = None):
        self.name = name
        self.values = channel_values
        self.brightness = brightness
        if all_bonds == None:
            self.all_bonds = []
        else:
            self.all_bonds = all_bonds
        self.committed = False
        self.final_bond = None
    def add_bond(self, bond):
        if bond not in self.all_bonds:
            self.all_bonds.append(bond)
        else:
            pass
           # raise Exception("Fluor " + self.name + "already had" + bond.name)
    
    def get_available_bonds(self):
        availabel_bonds = []
        for i in self.all_bonds:
            if i.state == "available":
                availabel_bonds.append(i)
        return availabel_bonds
    
    def __str__(self):
        return self.name
    def __repr__(self):
        return self.name
    def __eq__(self, other):
        if isinstance(other, Fluor):
            return self.name == other.name
    def __hash__(self):
        return hash(self.name)
    
class Marker:
    global visited_nucleotide_list
    def __init__(self, name: str, abundance: float, all_bonds = None):
        self.name = name
        self.abundance = abundance
        self.all_bonds = all_bonds
        if all_bonds == None:
            self.all_bonds = []
        else:
            self.all_bonds = all_bonds
        self.committed = False
        self.final_bond = None
    def add_bond(self, bond):
        if bond not in self.all_bonds:
            self.all_bonds.append(bond)
        else:
            pass
            #raise Exception("Marker " + self.name + "already had" + bond.name)
  
    def get_available_bonds(self):
        availabel_bonds = []
        for i in self.all_bonds:
            if i.state == "available":
                availabel_bonds.append(i)
        return availabel_bonds
    def __str__(self):
        return self.name
    def __repr__(self):
        return self.name
    def __eq__(self, other):
        if isinstance(other, Marker):
            return self.name == other.name
    def __hash__(self):
        return hash(self.name)
    
class Bond:
    def __init__(self, fluor_end: Fluor, marker_end: Marker, state: str):
        self.fluor_end = fluor_end
        self.marker_end = marker_end
        self.strength = np.log(fluor_end.brightness) * np.log(marker_end.abundance) / 10 
        self.state = state 
        self.name = self.marker_end.name + "---" \
        + self.fluor_end.name 
        
        # Update fluor and marker bond list 
        self.fluor_end.add_bond(self)
        self.marker_end.add_bond(self)
        self.intensity = None
        self.dna = None
    def __str__(self):
        return self.name        
    def __repr__(self):
        return self.name
    def __eq__(self, other):
        if isinstance(other, Bond):
            return self.name == other.name
    def __hash__(self):
        return hash(self.name)
    
class Fluor_strand:
    def __init__(self, all_fluors = None, ci = 0):
        if all_fluors == None:
            self.all_fluors = []
        else:
            self.all_fluors = all_fluors
        self.ci = ci
        self.fluor_strand_dict = {}
    def add_fluor(self, fluor: Fluor):
        if fluor not in self.all_fluors:
            self.all_fluors.append(fluor)
        else:
            pass
            #raise Exception("Fluor strand already has " + fluor.name)
    def sort_strand(self, small_first: bool = True):
        if small_first:
            self.all_fluors.sort(key=lambda x: len(x.all_bonds))
        else:
            self.all_fluors.sort(key=lambda x: len(x.all_bonds), reverse = True)    

    # Get the non-commited fluor which has the lowest number of available bonds
    def get_min_fluor(self):
        current_min = 100000
        min_fluor = None
        for i in self.all_fluors:
            available_bonds_num = len(i.get_available_bonds())
            if available_bonds_num < current_min and not i.committed:
                current_min = available_bonds_num
                min_fluor = i
        return min_fluor
    def generate_fluor_dict(self):
        for i in self.all_fluors:
            self.fluor_strand_dict[i.name] = i
    def __str__(self):
        name_string = "---Fluor Strand--- \n"
        for i in self.all_fluors:
            name_string += str(i.name) + ", "
        return name_string 
    def __repr__(self):
        name_string = "---Fluor Strand--- \n"
        for i in self.all_fluors:
            name_string += str(i.name) + ", "
        return name_string      

class Marker_strand:
    def __init__(self, all_markers = None):
        if all_markers == None:
            self.all_markers = []
        else:
            self.all_markers = all_markers
        self.marker_strand_dict = {}
    #    for i in all_markers:
#            self.marker_strand_dict[i.name] = i
    def add_marker(self, marker: Marker):
        if marker not in self.all_markers:
            self.all_markers.append(marker)
        else:
            pass
            #raise Exception("Fluor strand already has " + fluor.name)
    def sort_strand(self, small_first: bool = True):
        if small_first:
            self.all_markers.sort(key=lambda x: len(x.all_bonds))
        else:
            self.all_markers.sort(key=lambda x: len(x.all_bonds), reverse = True)    
            
    # Get the non-commited marker which has the lowest number of available bonds
    def get_min_marker(self):
        current_min = 100000
        min_marker = None
        for i in self.all_markers:
            available_bonds_num = len(i.get_available_bonds())
            if available_bonds_num < current_min and not i.committed:
                current_min = available_bonds_num
                min_marker = i
        return min_marker
    def generate_marker_dict(self):
        for i in self.all_markers:
            self.marker_strand_dict[i.name] = i
    def __str__(self):
        name_string = "---Marker Strand--- \n"
        for i in self.all_markers:
            name_string += str(i.name) + ", "
        return name_string  
    def __repr__(self):
        name_string = "---Marker Strand--- \n"
        for i in self.all_markers:
            name_string += str(i.name) + ", "
        return name_string  
    
class DNA:
    def __init__(self, fluor_strand: Fluor_strand, marker_strand: Marker_strand):
        self.fluor_strand = fluor_strand
        self.marker_strand = marker_strand
        self.all_bonds = []
        self.bond_formed_num = 0
        for i in marker_strand.all_markers:
            self.all_bonds.append(i.all_bonds)
        self.score = None
        self.bond_history = []
        self.marker_history = []
        self.fluor_history = []
        self.abundance_list = []
        self.brightness_list = []
    def __str__(self):
        bond_string = "---Bonds--- "
        for i in self.all_bonds:
            bond_string += "\n["
            for k in i:
                bond_string += str(k) +"_"+ str(k.state) + ", "
            bond_string += "]"
        dna_string = "[[[---DNA---]]]\n" + str(self.fluor_strand) + "\n" + \
        str(self.marker_strand) + "\n" + bond_string + "\n\n"
        return dna_string
    def __repr__(self):
        bond_string = "---Bonds--- "
        for i in self.all_bonds:
            bond_string += "["
            for k in i:
                if k.state == "formed":
                    bond_string += str(k) +"_"+ str(k.state) + ", "
            bond_string += "]"
        dna_string = "[[[---DNA---]]]\n" + str(self.marker_strand) + "\n" + \
        str(self.fluor_strand) + "\n" + bond_string + "\n\n"
       
        return bond_string

class Frontier:
    def __init__(self, comby_list):
        self.comby_list = comby_list
        self.reserve = [] 
        
    def length(self):
        return len(self.comby_list)
    
    def get_panel_list(self):
        return self.comby_list
    
    def get_panel_list_length(self):
        length_list = []
        for c in comby_list:
            length_list.append(c.length())
        return length_list
    
    def add_panel(self, comby):
        self.comby_list.append(comby)
        
    def hasPanel(self, comby):
        for i in self.get_comby_list():
            if i.equal(comby):
                return True
        return False 
    
        
    def __repr__(self):
        myRep = "Frontier = [\n"
        for comby in self.get_panel_list():
            myRep += comby.__repr__() + "\n"
        
        myRep += "]"
        return myRep
    
class Sample:
    def __init__(self, ID, cell_num:int = 100000, 
                dna: DNA = None, cells:list = None,\
                channel_num:int = 54): 
        self.id = ID
        self.dna = deepcopy(dna)
        self.cell_num = cell_num
        if cells == None:
            self.cells = []
            for i in range(cell_num):
                cell = Cell(ID = i, dna = dna)
                cell.channel_num = channel_num
                self.cells.append(cell)
        else:
            self.cells = cells
        self.components = defaultdict(list)
            
    def __str__(self):
        return "Sample_{0}({1} cells)".\
                format(self.id, self.cell_num)
    
    def __repr__(self):
        return self.__str__()
    
class Cell:
    def __init__(self, ID: int, dna: DNA):
        self.id = ID
        self.dna = dna
        self.added_fluors = []
        self.intensities = [] #same order as added flours
        self.min_abundance, self.max_abundance, \
        self.min_brightness, self.max_brightness = self.get_min_max_values()
        self.channel_num = 0
    
    def add_fluor(self, fluor: Fluor, intensity_factor: float):
        if fluor in self.added_fluors:
            return None
        self.added_fluors.append(fluor)
        # TO DISCUSS
        if fluor.name == "AfHuLym":
            self.intensities.append(intensity_factor)
        else:
            abundance = fluor.final_bond.marker_end.abundance
            brightness = fluor.brightness
            geo_mean_brightness = np.sqrt(self.min_brightness*self.max_brightness)
            geo_mean_abundance = np.sqrt(self.min_abundance*self.max_abundance)
            intensity = abundance/geo_mean_abundance * brightness/geo_mean_brightness * intensity_factor
            self.intensities.append(intensity)
        
    def emit_light(self, noise: float):
        emit_array = np.zeros(self.channel_num)
        for fluor, intensity in tuple(zip(self.added_fluors, self.intensities)):
            emit_array += fluor.values*intensity 
        #add noise #TO DISCUSS
        zero_channels = (emit_array == 0) * 1
        noise = np.random.normal(0, abs(noise*(emit_array)))
        return emit_array+noise
    
    def get_min_max_values(self):
        min_abundance, max_abundance = min(self.dna.abundance_list), max(self.dna.abundance_list)
        min_brightness, max_brightness = min(self.dna.brightness_list), max(self.dna.brightness_list)
        return min_abundance, max_abundance, min_brightness, max_brightness
    
    def __str__(self):
        cell_string = "Cell_"  + str(self.id) + "\n"
        for fluor, intensity in tuple(zip(self.added_fluors, self.intensities)):
            cell_string += str(fluor) + ": " \
                            + str(round(intensity)) \
                                + "\n"
        return cell_string 
    
    def __repr__(self):
        return self.__str__()
    
