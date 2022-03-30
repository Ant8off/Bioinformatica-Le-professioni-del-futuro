#=========================================================================================
#
#           Letture dei file Fasta
#
#=========================================================================================

#---------------------GENITORE 1

file1 = open('Fasta\genitore1\genitore1_alleleA.fa')                            #apertura del file
genitore1A = []

for linea in file1:                                             #ciclo per il riempimento della lista
    if not linea.startswith(">"):                               #se la riga non inizia per > la inserisce nella lista  
        genitore1A.append(linea.strip())

gen1A = "".join(genitore1A)                                     #conversione da lista a stringa

file2 = open('Fasta\genitore1\genitore1_alleleB.fa')                            #apertura file 2
genitore1B = []

for linea in file2:                                             #ciclo per il riempimento della lista        
    if not linea.startswith(">"):                               #se la riga non inizia per > la inserisce nella lista     
        genitore1B.append(linea.strip())

gen1B = "".join(genitore1B)                                     #conversione da lista a stringa
file1.close()
file2.close()
#----------------------GENITORE 2----------------

file3 = open('Fasta\genitore2\genitore2_alleleA.fa')                            #apertura del file
genitore2A = []

for linea in file3:                                             #ciclo per il riempimento della lista
    if not linea.startswith(">"):                               #se la riga non inizia per > la inserisce nella lista
        genitore2A.append(linea.strip())

gen2A = "".join(genitore2A)                                     #conversione da lista a stringa

file4 = open('Fasta\genitore2\genitore2_alleleB.fa')                            #apertura del file
genitore2B = []

for linea in file4:                                             #ciclo per il riempimento della lista
    if not linea.startswith(">"):                               #se la riga non inizia per > la inserisce nella lista
        genitore2B.append(linea.strip())

gen2B = "".join(genitore2B)                                     #conversione da lista a stringa
file3.close()
file4.close()
#--------------------------FIGLIO---------------------------

file6 = open('Fasta\Figlio\Figlio_alleleA.fa')                  #apertura del file
LFIA = []

for linea in file6:                                             #ciclo per il riempimento della lista
    if not linea.startswith(">"):                               #se la riga non inizia per > la inserisce nella lista
        LFIA.append(linea.strip())

FIA = "".join(LFIA)                                             #conversione da lista a stringa

file7 = open('Fasta\Figlio\Figlio_alleleB.fa')                               #apertura del file                              
LFIB = []

for linea in file7:                                             #ciclo per il riempimento della lista
    if not linea.startswith(">"):                               #se la riga non inizia per > la inserisce nella lista
        LFIB.append(linea.strip())

FIB = "".join(LFIB)                                             #conversione da lista a stringa
#print(len(FIB))

file6.close()
file7.close()
#========================================================================================
#
#               Controllo per le basi codificate dei genitori
#
#========================================================================================


#--------------------------GENITORE 1--------------------------

vcf = open('VCF\chr9_varianti.vcf')

split_lineGEN1 = []
varianti_codificateGEN1 = []
posizioniGEN1 = []
contatoreGEN1 = 0

for line in vcf:                                                                                    #Lettura delle righe nel file vcf
    if not line.startswith("#"):                                                                    #se la riga non inizia per # prende la base, la aggiunge ad una stringa e la rende in maiuscolo
        split_lineGEN1 = line.strip().split("\t")
        baseA = gen1A[int(split_lineGEN1[1])-1].upper()
        baseB = gen1B[int(split_lineGEN1[1])-1].upper()
    if baseA == split_lineGEN1[3]:                                                                  #controllo se le basi sono uguali o diverse   0 == |  1!= | . mancante
        baseA_codificata = "0"
    elif baseA==split_lineGEN1[4]:
        baseA_codificata = "1"
    else:
        baseA_codificata = "."
    if baseB == split_lineGEN1[3]:
        baseB_codificata = "0"
    elif baseB == split_lineGEN1[4]:
        baseB_codificata = "1"
    else:
        baseB_codificata = "."
    varianti_codificateGEN1.append(baseA_codificata + "/" + baseB_codificata)                       #aggiunge alla lista se == o !=    formato: "   0 / 1"
    if baseA_codificata == "1" or baseB_codificata == "1":                                          #aggiunge alla lista le varianti (almeno 1 1/1)
        posizioniGEN1.append(contatoreGEN1)
    contatoreGEN1 += 1
vcf.close()

#-------------------------------GENITORE 2------------------------------------------------------------------------------------

vcf = open('VCF\chr9_varianti.vcf')

split_lineGEN2 = []
varianti_codificateGEN2 = []
posizioniGEN2 = []
contatoreGEN2 = 0

for line in vcf:
    if not line.startswith("#"):                                                                    #Lettura delle righe nel file vcf
        split_lineGEN2 = line.strip().split("\t")                                                   #se la riga non inizia per # prende la base, la aggiunge ad una stringa e la rende in maiuscolo
        baseA = gen2A[int(split_lineGEN2[1])-1].upper()
        baseB = gen2B[int(split_lineGEN2[1])-1].upper()
    if baseA == split_lineGEN2[3]:                                                                  #controllo se le basi sono uguali o diverse   0 == |  1!= | . mancante
        baseA_codificata = "0"
    elif baseA==split_lineGEN2[4]:
        baseA_codificata = "1"
    else:
        baseA_codificata = "."
    if baseB == split_lineGEN2[3]:
        baseB_codificata = "0"
    elif baseB == split_lineGEN2[4]:
        baseB_codificata = "1"
    else:
        baseB_codificata = "."
    varianti_codificateGEN2.append(baseA_codificata + "/" + baseB_codificata)                       #aggiunge alla lista se == o !=    formato: "   0 / 1"
    if baseA_codificata == "1" or baseB_codificata == "1":                                          #aggiunge alla lista le varianti (almeno 1 1/1)
        posizioniGEN2.append(contatoreGEN2)
    contatoreGEN2 += 1
# print(varianti_codificateGEN2)
# print(posizioniGEN2)
vcf.close()

#-----------------------------------------Figlio-------------------------------------------------------------------

vcf = open('VCF\chr9_varianti.vcf')

split_lineFI = []
varianti_codificateFI = []
posizioniFI = []
contatoreFI = 0

for line in vcf:                                                                                    #Lettura delle righe nel file vcf
    if not line.startswith("#"):                                                                    #se la riga non inizia per # prende la base, la aggiunge ad una stringa e la rende in maiuscolo
        split_lineFI = line.strip().split("\t")
        baseA = FIA[int(split_lineFI[1])-1].upper()
        baseB = FIB[int(split_lineFI[1])-1].upper()
    if baseA == split_lineFI[3]:                                                                    #controllo se le basi sono uguali o diverse   0 == |  1!= | . mancante
        baseA_codificata = "0"
    elif baseA==split_lineFI[4]:
        baseA_codificata = "1"
    else:
        baseA_codificata = "."
    if baseB == split_lineFI[3]:
        baseB_codificata = "0"
    elif baseB == split_lineFI[4]:
        baseB_codificata = "1"
    else:
        baseB_codificata = "."
    varianti_codificateFI.append(baseA_codificata + "/" + baseB_codificata)                         #aggiunge alla lista se == o !=    formato: "   0 / 1"
    if baseA_codificata == "1" and baseB_codificata == "1":                                         #aggiunge alla lista le varianti (entrambi 1/1)
        posizioniFI.append(contatoreFI)
    contatoreFI += 1
vcf.close()
# print(varianti_codificateFI)
# print(posizioniFI)

#================================================================
#
# crea la lista dei "vcf"
#
#================================================================

vcf = open('VCF\chr9_varianti.vcf')

listapos = []
contpos = 0
contI = 0
listavcf = []

for line in vcf:
    listavcf.append(line)                           #aggiunge le linee del file vcf in una lista

vcf.close()

# print(listapos)
vcf = open('VCF\chr9_varianti.vcf')

#=======================================================
#
#aggiunge alla lista la posizione della base modificata
#
#=======================================================

o = 0
for i in range(0, len(listavcf)):                           #scorre la lunghezza della lista
    if o < len(posizioniFI):                                #se il contatore è minore della lunghezza
        if i == posizioniFI[o]:                             #se la linea == posizione della variazione 1/1 del figlio
            listapos.append(listavcf[i].split("\t"))        #aggiunge alla lista la riga dove è presente la variazione
            o += 1                                          #aumenta il contatore


lista = []
for i in range(0, len(listapos)):                           #scorre la lunghezza della lista delle posizioni delle variazioni
    lista.append(listapos[i][1])                            #aggiunge alla lista la posizione della variazione

#=======================================================
#
#Controllo se nelle righe file "bed" è presente "exon"
#
#=======================================================

gencode = open("Bed\gencode_chr9.bed")                      

exongen = []

for line in gencode:                                        #scorre il file bed e controlla che nelle righe sia presente la parola "exon"
    if "exon" in line:
        exongen.append(line.split("\t"))                    #aggiunge la linea in una lista


verifica = []

#==========================================================
#
#Controllo le basi modificate nel range presente nelle righe degli esoni
#
#=============================================================


for variante in lista:                                     #controlla che le posizioni salvate in listapos siano nel range presente nel file "bed" 
    variante = int(variante)     
    controllo = ""                         
    for line in exongen:
        
        start = int(line[1])
        end = int(line[2])

        if variante > start and variante < end and variante != controllo:   #stampa la posizione    
            controllo = variante 
            print("è presente una variante nella posizione:", variante)

