# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'genotyping.ui'
#
# Created by: PyQt5 UI code generator 5.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog 
from PyQt5.QtWidgets import QInputDialog

import sys
from matplotlib.patches import Patch
from matplotlib import colors as colorcode
from io import StringIO
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
import datetime
import random as rd
from Bio import SeqIO
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


class Ui_Form(object):
    def setupUi(self, Form,installationDirectory):
        Form.setObjectName("Form")
        Form.resize(794, 511)
        self.label = QtWidgets.QLabel(Form)
        self.label.setGeometry(QtCore.QRect(10, 10, 211, 20))
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(Form)
        self.label_2.setGeometry(QtCore.QRect(10, 90, 211, 20))
        self.label_2.setObjectName("label_2")
        self.label_3 = QtWidgets.QLabel(Form)
        self.label_3.setGeometry(QtCore.QRect(420, 10, 211, 20))
        self.label_3.setObjectName("label_3")
        self.label_4 = QtWidgets.QLabel(Form)
        self.label_4.setGeometry(QtCore.QRect(420, 90, 211, 20))
        self.label_4.setObjectName("label_4")
        self.label_5 = QtWidgets.QLabel(Form)
        self.label_5.setGeometry(QtCore.QRect(620, 90, 211, 20))
        self.label_5.setObjectName("label_5")
        self.label_6 = QtWidgets.QLabel(Form)
        self.label_6.setGeometry(QtCore.QRect(10, 190, 211, 20))
        self.label_6.setObjectName("label_6")
        self.inputFileEntry = QtWidgets.QLineEdit(Form)
        self.inputFileEntry.setGeometry(QtCore.QRect(10, 30, 231, 32))
        self.inputFileEntry.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.inputFileEntry.setObjectName("inputFileEntry")
        self.outputFolderEntry = QtWidgets.QLineEdit(Form)
        self.outputFolderEntry.setGeometry(QtCore.QRect(10, 110, 231, 32))
        self.outputFolderEntry.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.outputFolderEntry.setObjectName("outputFolderEntry")
        self.inputFileButton = QtWidgets.QPushButton(Form)
        self.inputFileButton.setGeometry(QtCore.QRect(250, 30, 112, 32))
        self.inputFileButton.setObjectName("inputFileButton")
        self.outputFolderButton = QtWidgets.QPushButton(Form)
        self.outputFolderButton.setGeometry(QtCore.QRect(250, 110, 112, 32))
        self.outputFolderButton.setObjectName("outputFolderButton")
        self.dbEntry = QtWidgets.QLineEdit(Form)
        self.dbEntry.setGeometry(QtCore.QRect(420, 30, 231, 32))
        self.dbEntry.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.dbEntry.setObjectName("dbEntry")
        self.dbButton = QtWidgets.QPushButton(Form)
        self.dbButton.setGeometry(QtCore.QRect(660, 30, 112, 32))
        self.dbButton.setObjectName("dbButton")
        self.detectionThresholdEntry = QtWidgets.QLineEdit(Form)
        self.detectionThresholdEntry.setGeometry(QtCore.QRect(420, 110, 161, 32))
        self.detectionThresholdEntry.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.detectionThresholdEntry.setObjectName("detectionThresholdEntry")
        self.numThreadsEntry = QtWidgets.QLineEdit(Form)
        self.numThreadsEntry.setGeometry(QtCore.QRect(620, 110, 151, 32))
        self.numThreadsEntry.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.numThreadsEntry.setObjectName("numThreadsEntry")
        self.logArea = QtWidgets.QTextEdit(Form)
        self.logArea.setGeometry(QtCore.QRect(10, 210, 541, 281))
        self.logArea.setObjectName("logArea")
        self.label_7 = QtWidgets.QLabel(Form)
        self.label_7.setGeometry(QtCore.QRect(580, 210, 201, 211))
        self.label_7.setText("")
        self.label_7.setPixmap(QtGui.QPixmap(installationDirectory+"GRACy_easyinstall/src/GUI/IconsFinal/Classification.jpg"))
        self.label_7.setObjectName("label_7")
        self.genotypeButton = QtWidgets.QPushButton(Form)
        self.genotypeButton.setGeometry(QtCore.QRect(580, 420, 201, 32))
        self.genotypeButton.setObjectName("genotypeButton")
        self.plotButton = QtWidgets.QPushButton(Form)
        self.plotButton.setGeometry(QtCore.QRect(580, 460, 201, 32))
        self.plotButton.setObjectName("plotButton")

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

         #Addition to GUI generated by Qt designer
        self.inputFileButton.clicked.connect(self.selectInputFile)
        self.outputFolderButton.clicked.connect(self.selectOutputFolder)
        self.plotButton.clicked.connect(self.plotGenotypes)
        self.genotypeButton.clicked.connect(self.runGenotyping)
        font = QtGui.QFont()
		font.setPointSize(11)
		Form.setFont(font)

    def selectInputFile(self):
        filename, __ = QFileDialog.getOpenFileName(None,"Select input file","./")
        self.inputFileEntry.setText(filename)


    def selectOutputFolder(self):
        folderName = QFileDialog.getExistingDirectory(None, "Select output folder","./")
        self.outputFolderEntry.setText(folderName)

    def runGenotyping(self):
        #Initialize variables
        orderedHyperLoci = ["rl5a","rl6","rl12","rl13","ul1","ul9","ul11","ul20","ul73","ul74","ul120","ul139","ul146"]

        #numReads = int(self.numReadsEntry.get())
        dbfile = installationDirectory+"src/scripts/genotyping/kmerDB/"+self.dbEntry.text()
        NumThreads = int(self.numThreadsEntry.text())
        inputFile = self.inputFileEntry.text()

        if not inputFile == "Please select file...." :
            if not self.outputFolderEntry.text() =="Please select a folder....":
                os.system("mkdir "+self.outputFolderEntry.text())
                datasets = []
                infile=open(inputFile)
                while True:
                    read1 = infile.readline().rstrip()
                    if not read1:
                        break
                    read2 = infile.readline().rstrip()
                    if not read2:
                        break

                    datasets.append((read1,read2))

                numDatasets = len(datasets)
                progressBarIncrement = 100/(numDatasets*20)
                step = 0
                for dataset in datasets:
                    read1 = dataset[0]
                    read2 = dataset[1]

                    #Extract path and files info
                    inputPath =  '/'.join(read1.split('/')[0:-1])
                    fileRoot1 = ((read1.split("/"))[-1])
                    if fileRoot1[-2:] == "fq":
                        fileRoot1 = fileRoot1[:-3]
                    else:
                        fileRoot1 = fileRoot1[:-6]

                    fileRoot2 = ((read2.split("/"))[-1])
                    if fileRoot2[-2:] == "fq":
                        fileRoot2 = fileRoot2[:-3]
                    else:
                        fileRoot2 = fileRoot2[:-6]

                    sampleRoot = fileRoot1[:-2]
                    outfile = open(fileRoot1+"_IDCard.txt","w")
                    logFile = open(fileRoot1+"_logFile.txt","w")
                    now = datetime.datetime.now()
                    logFile.write("Genotyping sample "+sampleRoot+"  started at "+now.strftime("%H:%M")+"\n")

                    fileList = "list"+str(rd.randint(0,1000000))+".txt"
                    listFile = open(fileList,"w")
                    listFile.write(read1+"\n"+read2+"\n")
                    listFile.close()

                    #Perform deduplication
                    now = datetime.datetime.now()
                    logFile.write("Deduplication for dataset "+fileRoot1+"  started at "+now.strftime("%H:%M")+"\n")

                    
                    
                    self.logArea.append("Genotyping dataset "+fileRoot1+"/"+fileRoot2+"\n")
                    self.logArea.repaint()
                    
                    
                    
                    


                    
                    self.logArea.append("Performing deduplication....")
                    self.logArea.repaint()
                    
                    
                    

                    #print "Performing deduplication...."
                    dedupFile1 = fileRoot1+"_noDup_1.fq"
                    dedupFile2 = fileRoot2+"_noDup_2.fq"
                    os.system(installationDirectory+"src/conda/bin/fastuniq -i "+fileList+" -t q -o "+dedupFile1+" -p "+dedupFile2)
                    os.system("rm -f "+fileList)

                    
                    
                    self.logArea.repaint()
                    
                    
                    

                    now = datetime.datetime.now()
                    logFile.write("Deduplication for dataset "+sampleRoot+"  ended at "+now.strftime("%H:%M")+"\n\n")
                    
                    
                    

                    #******************* Calculate average coverage for deduplicated reads   **************
                    now = datetime.datetime.now()
                    logFile.write("Reference coverage calculation for dataset "+sampleRoot+"  started at "+now.strftime("%H:%M")+"\n")

                    
                    self.logArea.append("Performing reference alignment....")
                    self.logArea.repaint()
                    
                    
                    
                    os.system(installationDirectory+"src/conda/bin/bowtie2 -1 "+dedupFile1+" -2 "+dedupFile2+" -x "+installationDirectory+"src/scripts/genotyping/fastaFiles/hcmvReference -S alignmenthsbfy43223.sam >null 2>&1")
                    
                    
                    self.logArea.repaint()
                    
                    
                    
                    
                    
                    

                    
                    self.logArea.append("Converting sam to bam....")
                    self.logArea.repaint()
                    
                    
                    
                    os.system(installationDirectory+"src/conda/bin/samtools view -bS -h alignmenthsbfy43223.sam > alignmenthsbfy43223.bam")
                    
                    
                    self.logArea.repaint()
                    
                    
                    
                    
                    
                    

                    
                    self.logArea.append("Sorting bam....")
                    self.logArea.repaint()
                    
                    
                    
                    os.system(installationDirectory+"src/conda/bin/samtools sort -o alignmenthsbfy43223_sorted.bam alignmenthsbfy43223.bam")
                    
                    
                    
                    
                    
                    
                    
                    
                    

                    
                    self.logArea.append("Calculating average coverage....")
                    self.logArea.repaint()
                    
                    
                    
                    os.system(installationDirectory+"src/conda/bin/samtools depth  alignmenthsbfy43223_sorted.bam  |  awk '{sum+=$3} END { print sum/NR}' >avCoverage.txt")
                    
                    
                    
                    
                    
                    
                    
                    

                    avCovFile = open("avCoverage.txt")
                    avCov = float(avCovFile.readline().rstrip())
                    avCovFile.close()
                    detectionTreshold = float(avCov*float(self.detectionThresholdEntry.text()))
                    
                    self.logArea.append("Average coverage for deduplicated reads: "+str(avCov)+"\n")
                    self.logArea.repaint()
                    
                    
                    
                    os.system("rm -f alignmenthsbfy43223.* avCoverage.txt "+fileList)
                    now = datetime.datetime.now()
                    logFile.write("Reference coverage calculation for dataset "+sampleRoot+"  ended at "+now.strftime("%H:%M")+"\n\n")

                    #Collect kmers from database
                    geneKmers = {}
                    kmerdbfile = open(dbfile)
                    kmerdbfile.readline()
                    while True:
                        line = kmerdbfile.readline().rstrip()
                        if not line:
                            break
                        fields = line.split("\t")
                        if not (fields[0],fields[1]) in geneKmers:
                            geneKmers[(fields[0],fields[1])] = []
                            kmerseqs = fields[2].split(",")
                            for item in kmerseqs:
                                if not len(item) == 0:
                                    kmerLength = len(item)
                                    geneKmers[(fields[0],fields[1])].append(item)


                    #Get sequences in memory
                    reads = []
                    numSeq = 0
                    overallGeneInfo = {}
                    
                    self.logArea.append("Caclculate kmer frequencies for file "+dedupFile1+"....")
                    self.logArea.repaint()
                    
                    
                    
                    os.system(installationDirectory+"src/conda/bin/jellyfish count -m 17 -s 100M -t 8 -C "+dedupFile1+" -o "+dedupFile1+"_kmerCount.jf")

                    
                    
                    
                    
                    
                    
                    
                    

                    numSeq = 0
                    
                    self.logArea.append("Caclculate kmer frequencies for file "+dedupFile2+"....")
                    self.logArea.repaint()
                    
                    
                    

                    os.system(installationDirectory+"src/conda/bin/jellyfish count -m 17 -s 100M -t 8 -C "+dedupFile2+" -o "+dedupFile2+"_kmerCount.jf")

                    
                    
                    
                    
                    


                    
                    self.logArea.append("Merging kmer files....")
                    self.logArea.repaint()
                    
                    
                    

                    os.system(installationDirectory+"src/conda/bin/jellyfish merge "+dedupFile1+"_kmerCount.jf "+dedupFile2+"_kmerCount.jf")
                    
                    
                    
                    
                    

                    
                    
                    

                    readsKmer = {}

                    
                    self.logArea.append("Loading kmers for sequences in "+dedupFile1+" into memory....")
                    self.logArea.repaint()
                    
                    
                    
                    for seq_record in SeqIO.parse(dedupFile1,"fastq"):
                        sequence = str(seq_record.seq)
                        revSequence = str(seq_record.seq.reverse_complement())
                        reads.append(str(seq_record.seq))
                        reads.append(str(seq_record.seq.reverse_complement()))
                        for a in range(0,len(sequence)-kmerLength+1):
                            if not sequence[a:a+kmerLength] in readsKmer:
                                readsKmer[sequence[a:a+kmerLength]] = []
                            readsKmer[sequence[a:a+kmerLength]].append(str(seq_record.id))

                        for a in range(0,len(revSequence)-kmerLength+1):
                            if not revSequence[a:a+kmerLength] in readsKmer:
                                readsKmer[revSequence[a:a+kmerLength]] = []
                            readsKmer[revSequence[a:a+kmerLength]].append(str(seq_record.id))

                        numSeq +=1
                        #if numSeq == int(numReads/2):
                        #    break

                    
                    
                    
                    
                    


                    
                    self.logArea.append("Loading kmers for sequences in "+dedupFile2+" into memory....")
                    self.logArea.repaint()
                    
                    
                    

                    for seq_record in SeqIO.parse(dedupFile2,"fastq"):
                        sequence = str(seq_record.seq)
                        revSequence = str(seq_record.seq.reverse_complement())
                        reads.append(str(seq_record.seq))
                        reads.append(str(seq_record.seq.reverse_complement()))
                        for a in range(0,len(sequence)-16):
                            if not sequence[a:a+17] in readsKmer:
                                readsKmer[sequence[a:a+17]] = []
                            readsKmer[sequence[a:a+17]].append(str(seq_record.id))

                        for a in range(0,len(revSequence)-16):
                            if not revSequence[a:a+17] in readsKmer:
                                readsKmer[revSequence[a:a+17]] = []
                            readsKmer[revSequence[a:a+17]].append(str(seq_record.id))

                        numSeq +=1
                        #if numSeq == int(numReads/2):
                        #    break

                    
                    
                    
                    
                    


                    now = datetime.datetime.now()
                    logFile.write("Genotyping for dataset "+sampleRoot+"  started at "+now.strftime("%H:%M")+"\n")

                    #Start searching for signatures
                    for gene in orderedHyperLoci:
                        
                        self.logArea.append("Genotyping gene "+gene+" for sample "+sampleRoot+"....")
                        self.logArea.repaint()
                        
                        
                        

                        matchedReads = {}
                        #Collect specific kmers for the genotypes of this gene
                        specificKmerGroup = {}
                        for item in geneKmers:
                            if item[0] == gene:
                                if not item[1] in specificKmerGroup:
                                    specificKmerGroup[item[1]] = []
                                for seqs in geneKmers[item]:
                                    specificKmerGroup[item[1]].append(seqs)




                        countSeq = {}
                        totCount = 0
                        numMatchedKmers = {}
                        for gr in specificKmerGroup:


                            if not gr in countSeq:
                                countSeq[gr] = 0

                            if not gr in matchedReads:
                                matchedReads[gr] = set()

                            if not gr in numMatchedKmers:
                                numMatchedKmers[gr] = 0

                            command = installationDirectory+"src/conda/bin/jellyfish query mer_counts_merged.jf "
                            for querySeq in specificKmerGroup[gr]:
                                command += querySeq
                                command += " "
                            command += " >counts.txt"

                            os.system(command)

                            countFile = open("counts.txt")


                            while True:
                                countLine = countFile.readline().rstrip()
                                if not countLine:
                                    break
                                countFields = countLine.split(" ")
                                if int(countFields[1])>0:
                                    numMatchedKmers[gr] +=1
                                    for item in readsKmer[countFields[0]]:
                                        matchedReads[gr].add(item)
                            countFile.close()
                            logFile.write(gene+"\t"+"For genotype "+gr+" there are "+str(len(matchedReads[gr]))+" reads matching\n")
                            if len(matchedReads[gr])>=detectionTreshold:
                                countSeq[gr] = len(matchedReads[gr])
                                totCount = totCount + len(matchedReads[gr])

                        
                        
                        
                        
                        

                        outfile.write(gene)


                        for gr in countSeq:
                            if countSeq[gr]>0:
                                percentage = float(countSeq[gr])/float(totCount)
                                outfile.write("\t"+gr+"\t"+str(percentage)+"\t"+str(totCount)+"\t"+str(numMatchedKmers[gr]))
                                
                                self.logArea.append("Found signature for genotype "+gr+" with a percentage of "+(str(percentage*100))[:5]+"%\n")
                                self.logArea.repaint()
                                
                                
                                

                        outfile.write("\n")
                        
                        
                        



                    now = datetime.datetime.now()
                    logFile.write("Genotyping sample "+fileRoot1+"  ended at "+now.strftime("%H:%M")+"\n")
                    logFile.close()
                    outfile.close()
                    os.system("rm -f "+dedupFile1+" "+dedupFile2+" *.jf")
                    os.system("rm -f alignmenthsbfy43223_sorted.bam null counts.txt")
                    os.system("mv *IDCard.txt "+self.outputFolderEntry.text())
                    os.system("mv *_logFile.txt "+self.outputFolderEntry.text())
            else:
                
                self.logArea.append("Please select a valid output folder to start genotyping....\n")
                self.logArea.repaint()
                
                
                

        else:
            
            self.logArea.append("Please select a valid input file to start genotyping....\n")
            self.logArea.repaint()
            
            
            
        
    def plotGenotypes(self):
        genes = ["ul120","ul9","rl12","rl13","ul1","ul11","ul139","rl5a","rl6","ul20","ul146","ul73","ul74"]
        colorDict = {"G1":"blue","G1A":"lightblue","G1B":"darkblue","G1C":"deepskyblue","G2":"red","G2A":"lightcoral","G2B":"darkred","G3":"salmon","G3A":"green","G3B":"darkgreen","G4":"lightgreen","G4A":"yellow","G4B":"darkseagreen","G4C":"pink","G4D":"maroon","G5":"magenta","G6":"orchid","G7":"purple","G8":"silver","G9":"blueviolet","G10":"darkcyan","G11":"navy","G12":"black","G13":"gray","G14":"lightgray"}
        
        if not self.outputFolderEntry.text() == "Please select a folder....":
            os.system("mkdir "+self.outputFolderEntry.text())
            packet = StringIO()
            numberOfPlots,__ = QFileDialog.getOpenFileNames(None,"Select IDcard file(s)","./")

            outputFileName, __ = QInputDialog.getText(None, 'Text Input Dialog', 'Enter plot name:')
 
            outputFileName = outputFileName+".png"

            
            self.logArea.append("Preparing plot....")


            shift = 0
            fig = plt.figure(figsize=(10,10))
            plot1 = fig.add_subplot(111)

            for sample in numberOfPlots:
                inputFile = sample
                hg = {}
                infile = open(inputFile)
                while True:
                    line = infile.readline().rstrip()
                    if not line:
                        break
                    fields = line.split("\t")
                    if not fields[0] in hg:
                        hg[fields[0]] = []
                    for a in range(1,len(fields),+4):
                        hg[fields[0]].append((fields[a],fields[a+1]))

                sizes = []
                colors = []

                for gene in genes:
                    for item in hg[gene]:
                        colors.append(colorDict[item[0]])
                        sizes.append(float(item[1]))
                    colors.append("white")
                    sizes.append(0.05)


                plot1.pie(sizes,radius=1-shift,colors=colors)
                plot1.pie([1],radius=0.95-shift,colors=["white"])
                
                shift = shift + 0.07

            yLabel = 0.5
            for sample in numberOfPlots:
                plt.text(-0.67,yLabel,(sample.split("/"))[-1],fontsize=14)
                yLabel = yLabel - 0.07

            legend_elements = []
            for item in colorDict:
                legend_elements.append(Patch(facecolor=colorcode.to_rgba(colorDict[item]),label=item))

            plt.legend(handles=legend_elements,loc=0,bbox_to_anchor=(1.12, 0.9))



            plot1.text(1,0.16,"UL120",rotation=14,fontweight='bold',fontsize=12,fontname="arial")
            plot1.text(0.80,0.61,"UL9",rotation=37,fontweight='bold',fontsize=12,fontname="arial")
            plot1.text(0.4,0.93,"RL12",rotation=62,fontweight='bold',fontsize=12,fontname="arial")
            plot1.text(-0.15,1.02,"RL13",rotation=-85,fontweight='bold',fontsize=12,fontname="arial")
            plot1.text(-0.66,0.84,"UL1",rotation=-60,fontweight='bold',fontsize=12,fontname="arial")
            plot1.text(-1.05,0.46,"UL11",rotation=-30,fontweight='bold',fontsize=12,fontname="arial")
            plot1.text(-1.22,-0.05,"UL139",rotation=3,fontweight='bold',fontsize=12,fontname="arial")
            plot1.text(-1.05,-0.57,"RL5A",rotation=24,fontweight='bold',fontsize=12,fontname="arial")
            plot1.text(-0.67,-0.94,"RL6",rotation=50,fontweight='bold',fontsize=12,fontname="arial")
            plot1.text(-0.17,-1.16,"UL20",rotation=75,fontweight='bold',fontsize=12,fontname="arial")
            plot1.text(0.35,-1.13,"UL146",rotation=-70,fontweight='bold',fontsize=12,fontname="arial")
            plot1.text(0.75,-0.81,"UL73",rotation=-45,fontweight='bold',fontsize=12,fontname="arial")
            plot1.text(0.97,-0.36,"UL74",rotation=-14,fontweight='bold',fontsize=12,fontname="arial")

            fig.savefig(self.outputFolderEntry.text()+"/"+outputFileName,dpi=300)

            

            self.logArea.append("Plot ready!\n")


            can = None
            new_pdf = None
            existing_pdf = None
            output = None
            page = None
            outputStream = None
        else:
            self.logArea.append("Please select a valid output folder to start genotyping....\n")
        


    

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Genotyping tool"))
        self.label.setText(_translate("Form", "Input file"))
        self.label_2.setText(_translate("Form", "Output folder"))
        self.label_3.setText(_translate("Form", "kmer database"))
        self.label_4.setText(_translate("Form", "Detection threshold"))
        self.label_5.setText(_translate("Form", "Number of threads"))
        self.label_6.setText(_translate("Form", "Log window"))
        self.inputFileEntry.setText(_translate("Form", "No file selected"))
        self.outputFolderEntry.setText(_translate("Form", "No folder selected"))
        self.inputFileButton.setText(_translate("Form", "Open file"))
        self.outputFolderButton.setText(_translate("Form", "Open folder"))
        self.dbEntry.setText(_translate("Form", "mainDB_seqs_filtered.txt"))
        self.dbButton.setText(_translate("Form", "Change DB"))
        self.detectionThresholdEntry.setText(_translate("Form", "0.02"))
        self.numThreadsEntry.setText(_translate("Form", "8"))
        self.genotypeButton.setText(_translate("Form", "Genotype!"))
        self.plotButton.setText(_translate("Form", "Plot"))

    

    

if __name__ == "__main__":
    import sys
    installationDirectory = sys.argv[1]
    app = QtWidgets.QApplication(sys.argv)
    Form = QtWidgets.QWidget()
    ui = Ui_Form()
    ui.setupUi(Form,installationDirectory)
    Form.show()
    sys.exit(app.exec_())
