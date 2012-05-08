#!/usr/bin/python
#-*-coding:utf-8-*-

'''
Created on 01/02/2011

@author: juliana
'''
import sys
import re
from urllib import urlopen


if __name__ == '__main__':

	try:
		archivo_original = sys.argv[1]
		destino = sys.argv[2]
		if destino[-1] != "/":
			destino = destino + "/"
	except IndexError:
		print('Falta el nombre del archivo o del destino: listar.py archivo.txt carpeta de destino')
		sys.exit(0)

	f = file( archivo_original , 'r')
	
	for nombreMotivo in f:
		if nombreMotivo != '\n': 
			nombreMotivo = nombreMotivo[:-1]
			elm_web = 'http://elm.eu.org/elmPages/' + nombreMotivo + '.html'
			pag = urlopen(elm_web).read()
			listaUID = re.findall( r'<a href="http://www.uniprot.org/uniprot/(\S+)">\1</a>', pag, re.IGNORECASE)
			nombre_archivo = destino + nombreMotivo + '_Ejemplos_UniprotID.txt'
			nuevo_archivo = open ( nombre_archivo , 'w')
			
			uniprotID = []
			
			for ej_elm in listaUID:
				if ej_elm not in uniprotID:
					uniprotID.append(ej_elm)
			
			for ej_elm in uniprotID:
				nuevo_archivo.write(ej_elm+"\n\n")
			nuevo_archivo.close()

	f.close()

