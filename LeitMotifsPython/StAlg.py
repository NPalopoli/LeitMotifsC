#Imports
from pylab import * #@UnusedWildImport
from scipy import * #@UnusedWildImport
from numpy import * #@UnusedWildImport
from math import * #@UnusedWildImport
import sets #@UnusedImport
from Bio import SeqIO #@UnresolvedImport
from weblogolib import *
from corebio import *
from scipy import stats
from copy import copy
import commands
import time

class Stormo(object):
	
	def _isOverlapping(self, p, M, l):
		for m in M:
			# Si pertenecen a la misma secuencia
			if (m[1] == p[1]):
				# Si overlappean
				if (m[0] < p[0] and m[0]+l >= p[0]) or (p[0] < m[0] and p[0]+l >= m[0]):
					return True
		return False		
	
	def _baseToPearsonVector(self, base):
	# El vector correspondiente a la base
		g = lambda x,b: 1.0 if x==b else 0.0
		return map(lambda x: g(x,base), self.bases)
	
	def _kwordToPearsonVector(self, kword):
	# Convierto una kword en un vector con unos en sus bases
		res = []
		for base in kword:
			for p in self._baseToPearsonVector(base):
				res += [p]
		return res
	
	def _resultsToPearson(self):
	# Convierto las palabras en vectores de apariciones de pearson	
		res = [[] for i in range(len(self.matrix))] #@UnusedVariable
		ind = 0
		for e in self.matrix:
			for e1 in e:
				if (res[ind] == []):
					# Si no puse nada solo agrego la matriz con unos en los lugares donde esten las bases
					res[ind] = self._kwordToPearsonVector(e1) 
				else:
					# Sino sumo contando apariciones de cada base en cada posicion
					res[ind] = map(sum, zip(res[ind],self._kwordToPearsonVector(e1)))
			ind += 1
		return res
	
	def _printM(self, M):
	# Imprime matrices de una manera mas amigable
		for e in M:
			print e
		print ""

	def _groupByPearson(self, pearsonfact):
		j = self._resultsToPearson()
		ind = 0
		# Empiezo de la posicion cero
		while ind < len(self.matrix):
			# Voy a comparar si desde estoy parado en adelante hay alguna que se parezca
			indcmp = ind+1
			# Hasta llegar al final
			while indcmp < len(self.matrix):
				# Si son lo suficientemente parecidas
				if (stats.pearsonr(j[ind], j[indcmp])[0] > pearsonfact):
					# Las junto
					self.matrix[ind] += self.matrix[indcmp] 
					self.positions[ind] += self.positions[indcmp]

					# Y borro la vieja		
					self.matrix.pop(indcmp)
					self.scores.pop(indcmp)
					self.scoreHist.pop(indcmp)
					self.positions.pop(indcmp)
					self._deleteRepeated()

					# Y recalculo las frecuencias
					j = self._resultsToPearson()

				else:		  
					indcmp+=1
			ind+=1
			
	def _deleteRepeated(self):
	# Ahora que las junte voy a eliminar las apariciones repetidas
		i2 = 0
		# Por cada una de las secuencias
		for m in self.positions:
			# Comienzo desde el final
			i = len(m) - 1
			# Mientras tenga algo por recorrer
			while (i>0 and i < len(self.matrix[i2])):
				# Si hay repetidos de lo que estoy buscando
				if (m.count(m[i]) > 1):
					# Si la posicion en la que estoy parado no es la primera aparicion
					if (m.index(m[i]) != i):
						# Borro la posicion en la que estoy parado
						m.pop(i)
						self.matrix[i2].pop(i)
				i-=1
			# Aumento uno en matrix
			i2+=1

	def _UploadSeqsFromFile(self, file):
	# Para levantar desde archivo FASTA
		f = open(file)
		for seq_record in  SeqIO.parse(f, "fasta"):
			self.seqs += [str(seq_record.seq)]
			self.seqnames += [seq_record.description]
		   
	def _SeqToMat(self):
	#Secuencia -> Matriz
		res = []
		self.totallength = 0
		for s in self.seqs:
			res += [list(s)]
			self.totallength += len(list(s)) 
		return res
	
	def _freqM(self,M, B):
	# Frecuencia matricial
		res = []
		for vec in M:
			res += [self._freq(vec,B)]
		return res
 
	def _freq(self, V, B):
	# Frecuencia de cada base en un vector
		res = []
		i = 0
		for base in B:
			res += [round(float(V.count(base)) / float(len(V)),2)]
			i+=1
		return res
 
	def _infM(self, M, B):
	# Informacion matricial
		res = []
		for vec in M:
			res += [self._inf(vec,B)]
		return res

	def _inf(self, F, P):
	# Informacion dadas frecuencias y probabilidades
		res = [4.32]
		for i in range(len(F)):
			if (F[i] != 0):
				res += [F[i]*log(F[i],2)]
		return res
	
	def _sumM(self, M, l):
	# Suma de vectores de una matriz
		res = []
		for v in M:
			res += [sum(v)]
		return map(lambda x: x - 19 / (2 * log(2) * l),res) 

	def _sum(self, V):
	# suma de un vector
		res = 0
		for i in range(len(V)):
				res += V[i]	 
		return res
	
	def _score(self, V):
		return sum(V)
	
	def _kwords(self, S,k):
	#Crea las matrices iniciales
		res = []
		for i in range(len(S) - k + 1):
			kword= S[i:k+i]
			res = res + [kword]
		return res
	
	def _preCalcule(self, M):
		# Creo un diccionario vacio
		self.precalculed = {}
		# Transpongo la matriz
		transposed = (transpose(M)).tolist()
		i = 0
	 	for e in transposed:
			i+=1
			for b in self.bases:
				n = e + [b]
				n = (transpose(n)).tolist()
				n = [[elem] for elem in n]
				self.precalculed[(b,i)] = self._getScore(n)

	def _getPreScore(self,kword):
	#Obtiene el score precalculado
		score = 0
		for i in range(len(kword)):
			# Al score le sumo el precalculo 
			# para la posicion y para las bases que intento agregar
			score += self.precalculed[kword[i], i+1]
	  	return score
   
	def _getScore(self, M):
	#Obtiene el score de informacion de una matriz
		transposed = (transpose(M)).tolist()
		frequency = self._freqM(transposed, self.bases)
		information = self._infM(frequency, self.freq_totales)
		summatory = self._sumM(information, len(M))
		return self._score(summatory)
	
	def _filterCorrelations(self):
		i = 0
		frequency = []
		for m in self.matrix:
			transposed = (transpose(m)).tolist()
			frequency += [self._freqM(transposed, self.bases)]
		return frequency
			
	def _shortestToTop(self, M):
		min = len(M[0])
		index = 0
		actual = 0
		for m in M:
			if len(m) < min:
				min = len(m)
				index = actual
			actual += 1
		aux = M[0]
		auxn =self.seqnames[0]
		M[0] = M[index]
		self.seqnames[0]=self.seqnames[index]
		M[index] = aux
		self.seqnames[index] = auxn
		return M

	def _sortByScores(self):
		for n in range(len(self.scores)): #@UnusedVariable
			temp = 0
			for i in range(1, len(self.scores)): 
				if self.scores[i] > self.scores[i-1]:
					temp = self.scores[i]
					temps = self.scoreHist[i]
					tempm = self.matrix[i]
					tempp = self.positions[i]
					
					self.scores[i] = self.scores[i-1]
					self.scoreHist[i] = self.scoreHist[i-1]
					self.matrix[i] = self.matrix[i-1]
					self.positions[i] = self.positions[i-1]
					
					self.scores[i-1] = temp 
					self.scoreHist[i-1] = temps 
					self.matrix[i-1] = tempm
					self.positions[i-1] = tempp
					
	def __init__(self, name, file):
	# Constructor
		self.bases = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
		self.freq_totales = [(1.0/20.0) for i in range (20)]
		self.name = name
		self.seqs = []
		self.seqnames = []
		self.positions = []
		self.scores = []
		self._UploadSeqsFromFile(file)
		self.matrix = []
	
	def _kwordsFiltered(self, Seqs, k):
		#Crea las matrices iniciales
		res = []
		for S in Seqs:
			if len(S) > k:
				kword= map(lambda e:e, S[0:k])
				for i in range(len(S) - k + 2):
					res = res + [kword]
					kword = map(lambda e:e, S[i:k+i])
		return res
	
	def runExhaustive(self, length, plotrange, protein):
	# Algoritmo exhaustivo
	
		# Paso de secuencia a matriz
		dataset = self._SeqToMat()
		lengthIni = self.totallength
		print "lengthIni: " + str(lengthIni)
		# Cargo todas las kwords y las aplano
		allkwords = reduce(lambda x, y: x+y, [self._kwords(ds, length) for ds in dataset])
		# Armo la matriz de posiciones
		allpositions = [(e+1,dataset.index(ds))  for ds in dataset for e in range(len(self._kwords(ds, length)))]
		# Las matrices iniciales son las primeras x kwords
		initialMatrices = allkwords
		# Defino la ubicacion de las posiciones iniciales
		self.positions = map(lambda x: [x], allpositions)
		# Armo una matriz a partir de una secuencia
		self.matrix = map(lambda x: [x], initialMatrices)
		# Seteo en blanco los scores
		self.scores = [self._getScore(self.matrix[0]) for e in range(len(initialMatrices))]
		# Historia de informacion
		self.scoreHist = [[] for e in range(len(initialMatrices))]
		j=0
		print "Preprocessed :", time.strftime("%H:%M:%S", time.localtime()),time.strftime("%d%b%Y", time.localtime())
		print ""
		print ""
		while j <= len(allkwords):
			# Agrego al analisis el siguiente conjunto de kwords
			newsequence = allkwords
			i=0
			# A cada paso agregare todas las kwords pertenecientes al conjunto inicial
			newmatrix = copy(self.matrix)
			for m in newmatrix:
				info = -1000.0
				maxpos = 0
				max = [] 
				pos = 0
				prevScore = self._getScore(m)
				# Aqui no busco el maximo score, sino que agrego todas al analisis
				for n in newsequence:
					# Tengo que cuidarme que no exista esa kword en la alineacion
					if self.positions[i].count(allpositions[pos]) == 0:
						# Secuencia que agregare al analisis
						actual = m + [n]
						# Score de esa secuencia
						score = self._getScore(actual) 
						# Agrego la secuencia al analisis
						self.matrix += [actual]
						self.scores += [score]
						self.positions += [self.positions[i] + [allpositions[pos]]]
					pos+=1
					#Esto es solo para imprimir el avance
				i+=1
			j+=1
			# Me guardo paso a paso cual es la informacion
			for i in range(len(self.scores)): 
				self.scoreHist += [self.scores[i] + self.scores[len(self.scores)-1]]
			print "Total : " + str(len(allkwords)+1) + " Actual progress: " + str(j)
			print "Step: " + str(j)
			print "Max score founded at :", time.strftime("%H:%M:%S", time.localtime()),time.strftime("%d%b%Y", time.localtime())
			print ""
			print ""
		# Ordeno scores
		self._sortByScores()
		i=0
		while i<1000:
			print self.matrix[i]
			print self.positions[i]
			print self.scores[i]
			print self.scoreHist[i]
			i = i+100
			raw_input("")
		# Escribo archivos
		self._writeFiles(protein, length, 0, plotrange, False, 0.0)
		
	def runLikeFullSeq(self, length, plotrange, protein, filter, pearson, pearsonfact, allow_overlapping):
	# Algoritmo principal
		# lenght: largo del patron buscado
		# lenghtIni: cantidad de kwords en las matrices iniciales
		# times: cantidad de veces que ofrezco cada secuencia
		
		# Paso de secuencia a matriz
		dataset = self._SeqToMat()
		# Cargo todas las kwords y las aplano
		allkwords = reduce(lambda x, y: x+y, [self._kwords(ds, length) for ds in dataset])
		# A cada paso contrastare todas los alineamientos contra todas las kwords del conjunto
		lengthIni = len(allkwords)
		# Armo la matriz de posiciones
		allpositions = [(e+1,dataset.index(ds))  for ds in dataset for e in range(len(self._kwords(ds, length)))]
		# Las matrices con alineamientos iniciales son las primeras lengthIni kwords
		initialMatrices = [allkwords[kw] for kw in range(0, lengthIni-1)]
		# Veo si en una ronda no cambian dejo de analizar, empiezo analizando todas
		changes = [True for kw in range(0, lengthIni)]
		# Defino la ubicacion de las posiciones iniciales
		self.positions = [[allpositions[kw]] for kw in range(0, lengthIni-1)]
		# Armo una matriz a partir de una secuencia
		self.matrix = map(lambda x: [x], initialMatrices)
		# Seteo en blanco los scores
		self.scores = [0 for e in range(len(initialMatrices))]
		# Historia de informacion
		self.scoreHist = [[] for e in range(len(initialMatrices))]
		j=0
		# Llegado este punto estoy listo para comenzar a correr el metodo
		print "Preprocessed :", time.strftime("%H:%M:%S", time.localtime()),time.strftime("%d%b%Y", time.localtime())
		print ""
		print ""
		# follow es el flag que me indica que esta ronda algun alineamiento cambio y debo seguir ejecutando
		follow = True
		while follow:
			# Agrego al analisis el siguiente conjunto de kwords
			newsequence = copy(allkwords)
			i = 0
			# Seteo follow en false y si alguna cambia paso el flag a true y eso marca que debo seguir
			follow = False
			# Comparo con cada una de las secuencias de la matriz de alineamientos
			for m in self.matrix:
				# Si este es uno de los alineamientos que han cambiado y que por lo tanto debo analizar
				if changes[i] == True:
					info = -1000.0
					pos = -1
					maxpos = 0
					max = [] 
					prevScore = self._getScore(m)
					# Precalculo scores
					self._preCalcule(m)
					# Buscando de todas las kwords que puedo alinear la que genera un nuevo alineamiento de maximo score
					for n in newsequence:
						# Tengo que cuidarme que no exista esa kword en la alineacion actualmente
						pos+=1
						if (self.positions[i].count(allpositions[pos])) == 0:
							if (not(allow_overlapping) and not(self._isOverlapping(allpositions[pos], self.positions[i], length))) or allow_overlapping:
								score = self._getPreScore(n)
								if score > info:
									max = m + [n]
									info = score
									maxpos = pos
					#Cuando la encuentro las seteo, si es que eso no hace que la informacion del alineamiento baje
					if prevScore * (1 - 0.05/(j+1)) < info:
						self.matrix[i] = max
						self.scores[i] = info
						self.positions[i] += [allpositions[maxpos]]
						#Si cambio alguna significa que tengo que seguir analizando este alineamiento en lo sucesivo
						follow = True
					else:
						#Si no agregue ya no agregare, asi que dejo de tenerla en cuenta
						changes[i] = False
				#Esto es solo para imprimir el avance
				i+=1
			j+=1
			# Me guardo paso a paso cual es la informacion de cada alineamiento para mantener un historial
			for i in range(len(self.scores)): 
				self.scoreHist[i] += [self.scores[i]]
			# Imprimo el avance
			print "Step: " + str(j)
			print "Max score founded at :", time.strftime("%H:%M:%S", time.localtime()),time.strftime("%d%b%Y", time.localtime())
			# Ordeno los resultados en base a la cantidad de informacion final
			print ""
			print ""
		# Filtro los resultados usando pearson si esa fue la eleccion del usuario
		self._sortByScores()
		print "Scores sorted at :", time.strftime("%H:%M:%S", time.localtime()),time.strftime("%d%b%Y", time.localtime())
		# Imprimo archivos y termine
		self._writeFiles(protein, length, j, plotrange, filter, 0.0)
		print "Files written at :", time.strftime("%H:%M:%S", time.localtime()),time.strftime("%d%b%Y", time.localtime())
		if pearson:
			self._groupByPearson(pearsonfact)
			print "Grouped by pearson at :", time.strftime("%H:%M:%S", time.localtime()),time.strftime("%d%b%Y", time.localtime())
			self._sortByScores()
			print "Scores sorted at :", time.strftime("%H:%M:%S", time.localtime()),time.strftime("%d%b%Y", time.localtime())
			# Imprimo archivos y termine
			self._writeFiles(protein, length, j, plotrange, filter, pearsonfact)
			print "Files written at :", time.strftime("%H:%M:%S", time.localtime()),time.strftime("%d%b%Y", time.localtime())
		
	def _plotAtPos(self, file, pos):
	# Dibuja el weblogo en el archivo dado de la posicion dada
		fout = open("weblogo/temp.fasta", 'w')
		secus=set()
		for e in range(len(self.positions[pos])):
			fout.write(">"+self.seqnames[self.positions[pos][e][1]]+"\n"+str(reduce(lambda x,y: x+y, self.matrix[pos][e])+"\n"))
			secus.add(self.seqnames[self.positions[pos][e][1]])
		fout.close()
		commands.getoutput('./weblogo/seqlogo -f weblogo/temp.fasta -o '+file+' -F PNG -c -e -Y -x Info_Total:'+str(round(self.scores[pos],2))+' -M -t "Kwords:'+str(len(self.matrix[pos]))+'-Secuencias:'+str(len(secus))+'" ')
			
	def _printPositions(self, file, pos):
	# Imprime para cada una de las secuencas a partir de que posicion comienza el patron
		fout = open(file, 'w')
		for e in range(len(self.positions[pos])):
			fout.write(str(self.positions[pos][e][0])+"\t"+reduce(lambda x,y: x+y, self.matrix[pos][e])+"\t"+self.seqnames[self.positions[pos][e][1]]+"\n")
			
	def _writeFiles(self, protein, length, time, plotrange, filter, pearsonfact):
		#Creacion de los archivos
		cla()
	        dir="/home/npalopoli/LeitMotifs/Resultados/"+str(pearsonfact)+"/"
#		dir="/home/npalopoli/LeitMotifs/Pruebas/"+str(pearsonfact)+"/"
		commands.getoutput("mkdir "+dir)
		for i in range(min(plotrange, len(self.matrix))):
			resdir=dir+protein+"-"+str(length)+"-"+str(time)+"-"+str(filter)+"/"
			commands.getoutput("mkdir "+resdir)
			# WebLogo
			self._plotAtPos(resdir+str(i), i)
			# Historial -> Archivo
			plot([self.scoreHist[i][j]  for j in range(len(self.scoreHist[i]))])
			savefig(resdir+str(i)+"Hist.png")
			# Creacion de un archivo con las posiciones donde esta el patron hallado
			self._printPositions(resdir+str(i)+".txt", i)

