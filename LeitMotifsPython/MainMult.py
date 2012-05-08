from StAlg import * #@UnusedWildImport
from VSL2 import * #@UnusedWildImport

# Defino la funcion que ejecutara las corridas
def Run(protein, length, plotrange, filterCDHIT, CDHIT_perc, filterVSL2, pearson, pearsonfact, exhaustive, allow_overlapping):
    print "Start :", time.strftime("%H:%M:%S", time.localtime()),time.strftime("%d%b%Y", time.localtime())

    # Parametros de la corrida
    work_file  = "/home/npalopoli/LeitMotifs/Filter4plus/Secuencias/TRG_ER_FFAT_1_Ejemplos_Secuencias.fasta"

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

Run("TRG_ER_FFAT_1_Ejemplos_Secuencias.fasta", 12, 100, True, 0.63, True, False, 0.98, False, False)
