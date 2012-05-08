#!/bin/bash
#======================================================================================
#
#         FILE: MainMult.sh
#
#        USAGE: MainMult.sh [-d] [-l] [-oD logfile] [-h] [starting directories]
#
#  DESCRIPTION: Run LeitMotifs for all FASTA alignments specified in file.
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Nicolas Palopoli
#      COMPANY: Universidad de Buenos Aires
#      CONTACT: nicopalo@gmail.com
#      VERSION: 1.0
#      CREATED: ---
#     REVISION: ---
#
#======================================================================================

if [ -z "$1" ]; then
	echo "Please specify path to file with list of FASTA files."
	exit
else
	if [ ! -e "$1" ]; then
		echo "Can't open list of FASTA files."
		exit
	else	
		LIST="$1"
	fi
fi

if [ ! -d "$2" ]; then
	echo "Please specify path to directory of FASTA files."
	exit
else
	DIR=`echo "$2" | sed 's/\/$//g'`
fi


MAIN1=$( cat <<MAINtmp
from StAlg import * #@UnusedWildImport
from VSL2 import * #@UnusedWildImport

# Defino la funcion que ejecutara las corridas
def Run(protein, length, plotrange, filterCDHIT, CDHIT_perc, filterVSL2, pearson, pearsonfact, exhaustive, allow_overlapping):
    print "Start :", time.strftime("%H:%M:%S", time.localtime()),time.strftime("%d%b%Y", time.localtime())

    # Parametros de la corrida
MAINtmp
)

MAIN2=$( cat <<MAINtmp

    #Me preparo para la eventual aplicacion de filtros
    t = Filter(work_file)

    if filterCDHIT:
        # Usando CD-Hit
        t.runCDHIT(CDHIT_perc)
        work_file = work_file + ".fil"
        t.filename = work_file
    if filterVSL2:
        # Usando el filtro VSL2
        t.runVSL2()
        work_file = work_file + ".fil"
        t.filename = work_file

    # Construyo el objeto con las secuencias
    s = Stormo(protein, work_file)

    # Dependiendo de la modalidad de corrida
    if exhaustive:
        s.runExhaustive(length, plotrange, protein)
    else:
        s.runLikeFullSeq(length, plotrange, protein, filterCDHIT, pearson, pearsonfact, allow_overlapping)

    # Imprimo la hora de finalizacion
    print "End :", time.strftime("%H:%M:%S", time.localtime()),time.strftime("%d%b%Y", time.localtime())

    return s

MAINtmp
)


while read line
do
	NAME=`echo "$line" | cut -d' ' -f 1`
	LEN=`echo "$line" | cut -d' ' -f 2`
        echo "$MAIN1" >MainMult.py
	echo "    work_file  = \"$DIR/$NAME\"" >>MainMult.py
	echo "$MAIN2" >>MainMult.py
	echo >>MainMult.py
	echo "Run(\"$NAME\", $LEN, 100, True, 0.63, True, False, 0.98, False, False)" >>MainMult.py
	chmod 777 MainMult.py
	python MainMult.py
done<"$LIST"

