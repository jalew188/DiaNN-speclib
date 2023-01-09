import struct
import pandas as pd
import numpy as np
import tqdm

def read_int32(f):
    return struct.unpack('i',f.read(4))[0]
def read_double(f):
    return struct.unpack('d',f.read(8))[0]
def read_float(f):
    return struct.unpack('f',f.read(4))[0]
def read_string(f):
    size = read_int32(f)
    if size == 0: return ""
    return f.read(size).decode('ascii')

def read_vector(f, type_size=4):
    size = read_int32(f)
    return f.read(size*type_size)

def read_objects(f, _type):
    size = read_int32(f)
    objects = [_type() for i in range(size)]
    for item in objects:
        item.read(f)
    return objects

def read_strings(f):
    size = read_int32(f)
    return [read_string(f) for i in range(size)]

class Protein:
    def __init__(self):
        self.swissprot = 1
        self.id = ""
        self.name = ""
        self.gene = ""
        self.name_idx = 0
        self.gene_idx = 0
    def read(self, f):
        self.swissprot = read_int32(f)
        size = read_int32(f)
        self.id = read_string(f)
        self.name = read_string(f)
        self.gene = read_string(f)
        self.name_idx = read_int32(f)
        self.gene_idx = read_int32(f)
        f.read(4*size)

class PG:
    def __init__(self):
        self.ids = ""
        self.names = ""
        self.genes = ""
    
    def read(self,f):
        size = read_int32(f)
        self.ids = read_string(f)
        self.names = read_string(f)
        self.genes = read_string(f)
        read_vector(f)
        read_vector(f)
        read_vector(f)
        f.read(4*size)

class DiannLibrary:
    def read_entries(self,f):
        size = read_int32(f)
        self.precursor_df = pd.DataFrame(
            columns = 'lib_index,charge,nAA,mz,irt,srt,unknown,im,iim'.split(',')
        )
        self.precursor_df['lib_index'] = [int(0)]*size
        self.precursor_df['charge'] = int(0)
        self.precursor_df['nAA'] = int(0)
        self.precursor_df['mz'] = 0.0
        self.precursor_df['irt'] = 0.0
        self.precursor_df['srt'] = 0.0
        self.precursor_df['unknown'] = 0.0
        self.precursor_df['im'] = 0.0
        self.precursor_df['iim'] = 0.0
        frag_start_idxes = [0]*size
        frag_stop_idxes = [0]*size
        precursor_names = [""]*size
        proteins = ['']*size
        genes = ['']*size
        proteotypics = [0]*size
        self.fragment_df = pd.DataFrame(
            columns = 'mz,intensity,charge,frag_type,index,loss'.split(',')
        )
        self.fragment_df['mz'] = self.fragment_df['mz'].astype(np.float64)
        self.fragment_df['intensity'] = self.fragment_df['intensity'].astype(np.float64)
        self.fragment_df['charge'] = self.fragment_df['charge'].astype(np.int8)
        self.fragment_df['frag_type'] = self.fragment_df['frag_type'].astype(np.int8)
        self.fragment_df['index'] = self.fragment_df['index'].astype(np.int8)
        self.fragment_df['loss'] = self.fragment_df['loss'].astype(np.int8)
        frag_info_list = []
        frag_start_idx = 0
        protein_ids = []
        for i in tqdm.tqdm(range(size)):
            self.precursor_df.loc[i,:] = struct.unpack('3i6f', f.read(36))
            frag_size = read_int32(f)
            frag_start_idxes[i] = frag_start_idx
            frag_stop_idxes[i] = frag_start_idx + frag_size
            frag_start_idx = frag_stop_idxes[i]
            for _ in range(frag_size):
                frag_info_list.append(struct.unpack('2f4b',f.read(12)))
            dc = read_int32(f)
            if dc > 0:
                struct.unpack('3i6f', f.read(36))
                frag_size = read_int32(f)
                for _ in range(frag_size):
                    struct.unpack('2f4b',f.read(12))
            entry_frag = read_int32(f)
            proteotypics[i] = read_int32(f)
            protein_id = read_int32(f)
            proteins[i] = self.PGs[protein_id].ids
            genes[i] = self.PGs[protein_id].genes
            precursor_names[i] = read_string(f)
            struct.unpack('3f',f.read(12)) # I don't know what are these entries ...
        self.precursor_df['frag_start_idx'] = frag_start_idxes
        self.precursor_df['frag_stop_idx'] = frag_stop_idxes
        self.precursor_df['precursor_name'] = precursor_names
        self.precursor_df['genes'] = genes
        self.precursor_df['proteins'] = proteins
        self.precursor_df['proteotypic'] = proteotypics
        self.fragment_df.loc[:] = frag_info_list
        
    def read(self,_io):
        if isinstance(_io, str):
            f = open(_io, 'rb')
        else:
            f = _io
        version = read_int32(f)
        if version > 0:
            gd = version
        else:
            gd = read_int32(f)

        gc = read_int32(f)
        ip = read_int32(f)

        name = read_string(f)
        self.fasta_names = read_string(f)

        self.proteins = read_objects(f, Protein)
        self.PGs = read_objects(f, PG)
        self.precursors = read_strings(f)
        self.names = read_strings(f)
        self.genes = read_strings(f)

        self.min_irt = read_double(f)
        self.max_irt = read_double(f)

        self.read_entries(f)
        if version <= -1:
            try:
                size = read_int32(f)
                f.read(4*size)
                print("Elution groups loaded")
            except Exception:
                print("Elution groups failed")
                pass
        if isinstance(_io, str):
            f.close()