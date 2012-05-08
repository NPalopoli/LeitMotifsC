
#Imports
from scipy import * #@UnusedWildImport
from numpy import * #@UnusedWildImport
from math import * #@UnusedWildImport
import matplotlib.pyplot as plt #@UnusedImport
import sets #@UnusedImport
from Bio import SeqIO #@UnresolvedImport
from weblogolib import *
from corebio import *
import commands

# Clase donde tendre metodos de filtrado que aplicare potencialmente sobre un archivo con secuencias
class Filter(object):
    
    # Armo la lista de secuencias presente en el archivo
    def _UploadSeqsFromFile(self, file):
    # Para levantar desde archivo FASTA
        for seq_record in  SeqIO.parse(file, "fasta"):
            self.seqs += [str(seq_record.seq)]
            #self.seqs = map(lambda element:reduce(lambda x,y: x+y, element), self.seqs)
            self.seqnames += [seq_record.description]
    
    # Armo un archivo en base a un conjunto de secuencias
    def _SeqToFile(self, seq, file):
    # Escribe una secuencia en archivo
        fout = open(file, "w")
        fout.write(seq)
        fout.close()
    
    
    # Carga las regiones desordenadas de una secuencia en una matriz
    def _LoadDisorderedRegions(self):
        # Me fijo cada caracter si esta o no ordenado
        fin = open("VSL2/out", "r")
        lines = fin.readlines()
        fin.close()
        lines = filter(lambda x: x.strip()[len(x.strip())-1] == '.' or x.strip()[len(x.strip())-1] == 'D', filter(lambda x: len(x.strip()) > 0, lines))
        lines = map(lambda x:(x[1],x[3]),map(lambda x: x.split(), lines))
        # Esto es para tomar los bordes de las regiones desordenadas
        for j in range(10): #@UnusedVariable
            for i in range((len(lines)-1)):
                if i > 0:
                    if lines[i-1][1] == "." and lines[i][1] == "D":
                        lines[i-1] = (lines[i-1][0], "D")
            i = len(lines)-2
            while i > 0:
                if i > 0:
                    if lines[i+1][1] == "." and lines[i][1] == "D":
                        lines[i+1] = (lines[i+1][0], "D")
                i-=1
                
        # Aplano y filtro las regiones desordenadas
        res = []
        actual = [lines[0][0]]
        anterior = lines[0][1]
        for l in lines:
            if l[1] == anterior:
                actual += l[0]
            else:
                if l[1] == '.':
                    res += [reduce(lambda x, y: x+y, actual)]
                actual = [l[0]]
            anterior = l[1]
        return res
        
    # Escribe el fasta, para cada secuencia pone todas sus regiones
    def _writeFasta(self):
        fout = open(self.filename, "w")
        for i in range(len(self.seqs)):
            for j in range(len(self.seqs[i])):
                fout.write(">"+self.seqnames[i]+"\n")
                fout.write(self.seqs[i][j]+"\n")
        fout.close()
        
    # Constructor de la clase
    def __init__(self, file):
        self.filename = file
        self.seqs = []
        self.seqnames = []
        self._UploadSeqsFromFile(file)
        
    # Arma el archivo filtrado con la herramienta VSL2
    def runVSL2(self):
        for i in range(len(self.seqs)):
            self._SeqToFile(self.seqs[i], "VSL2/temp")
            commands.getoutput("java -jar VSL2/VSL2.jar -s:VSL2/temp > VSL2/out")
            self.seqs[i] = self._LoadDisorderedRegions()
        self.seqs = map(lambda e: filter(lambda x: len(x) > 0, e), self.seqs)
        self.filename = self.filename+".fil"
        self._writeFasta()
        
    # Filtra el archivo utilizando cdhit
    def runCDHIT(self, porc):
        commands.getoutput("./cd-hit/cd-hit -i "+self.filename+" -o "+self.filename+".fil -c "+str(porc))
         
        
        
        
