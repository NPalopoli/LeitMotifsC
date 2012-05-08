from StAlg import * #@UnusedWildImport
from VSL2 import * #@UnusedWildImport

# Defino la funcion que ejecutara las corridas
def Run(protein, length, plotrange, filterCDHIT, CDHIT_perc, filterVSL2, pearson, pearsonfact, exhaustive, allow_overlapping):
    print "Start :", time.strftime("%H:%M:%S", time.localtime()),time.strftime("%d%b%Y", time.localtime())

    # Parametros de la corrida
    work_file  = "FASTA/"+protein+".fa"
#    work_file  = "in"+protein+".fasta"

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
    
# Ejecuto una corrida
#Run("AdHoc", 5, 100, False, 0.4, False, True, 0.98, False, True)
#Run("LIG_SH3_2", 6, 100, False, 0.63, False, True, 0.98, False, False)
#Run("LIG_14-3-3_1", 6, 100, False, 0.63, False, True, 0.98, False, False)
#Run("LIG_14-3-3_3", 6, 100, False, 0.63, False, True, 0.98, False, False)
# La linea a continuacion era activa en lo que me paso Leandro
#Run("LIG_AP_GAE_1", 7, 100, False, 0.63, False, False, 0.98, False, False)
#Run("LIG_Clathr_ClatBox_1", 5, 100, False, 0.63, False, True, 0.98, False, False)
#Run("LIG_CtBP", 5, 100, False, 0.63, False, True, 0.98, False, False)
#Run("LIG_EH1_1", 7, 100, False, 0.63, False, False, 0.98, False, False)
#Run("LIG_HP1_1", 7, 100, False, 0.63, False, True, 0.98, False, False)
#Run("LIG_NRBOX", 5, 100, False, 0.63, False, True, 0.98, False, False)
#Run("MOD_CMANNOS", 4, 100, False, 0.63, False, True, 0.98, False, False)
#Run("LIG_RGD", 4, 100, False, 0.63, False, True, 0.98, False, False)
#Run("01", 11, 100, False, 0.63, False, True, 0.98, False, False)
#Run("LxCxE", 5, 1000, True, 0.63, True, False, 0.98, False, False)
#Run("RB2", 5, 100, True, 0.63, True, False, 0.98, False, False)
#Run("RB", 5, 100, True, 0.63, True, False, 0.98, False, False)
#Run("LIG_TRAF2_1", 4, 100, True, 0.63, True, False, 0.98, False, False)
#Run("LIG_TRAF6", 3, 100, True, 0.63, True, False, 0.98, False, False)
#Run("LIG_PCNA", 8, 100, True, 0.63, True, False, 0.98, False, False)
#Run("LIG_Dynein_DLC8_1", 8, 100, True, 0.63, True, False, 0.98, False, False)

#Run("E7C", 5, 100, False, 0.63, False, False, 0.98, False, False)
#Run("E7C", 7, 100, False, 0.63, False, False, 0.98, False, False)
#Run("E7C", 9, 100, False, 0.63, False, False, 0.98, False, False)
#Run("E7C", 11, 100, False, 0.63, False, False, 0.98, False, False)

#Run("1N11", 5, 100, False, 0.63, False, False, 0.98, False, False)
#Run("1N11", 10, 100, False, 0.63, False, False, 0.98, False, False)
#Run("1N11", 20, 100, False, 0.63, False, False, 0.98, False, False)
#Run("1N11", 30, 100, False, 0.63, False, False, 0.98, False, False)
#Run("1N11", 40, 100, False, 0.63, False, False, 0.98, False, False)

#Run("Nico", 7, 100, False, 0.63, False, False, 0.98, False, False)

#Run("LIG_14-3-3_1", 6, 100, True, 0.63, True, True, 0.98, False, False)
#Run("LIG_Clathr_ClatBox_1", 5, 100, True, 0.63, True, True, 0.98, False, False)
#Run("LIG_Clathr_ClatBox_1", 5, 100, False, 0.63, False, True, 0.98, False, False)
#Run("Nico", 4, 100, False, 0.63, False, True, 0.98, False, False)
