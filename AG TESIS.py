# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 20:08:54 2022
@author: steve
"""

# Data:
Element=[['Link1','L1','L2','L3','L4','L5','L1Y','L2Y','L3Y','L4Y','L5Y']]
units='Kgf_cm_C'
Ke1Init=[[36162,72323,96431,114512,120538,36162,72323,96431,114512,120538]]   #kg/cm
PerObj=[0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02]      #Deriva Objetivo
ModalMassObj=[]
StructuralWeightObj=[]
StructureDemand_CapacityObj=[] 
ModelDirectory='C:\\Users\\steve\\Desktop\\CODE\\ETABS\\'
ProgramPath='C:\Program Files\Computers and Structures\ETABS 19\ETABS.exe'
ModelName0='MODELO_EVALUAR COL BRBF'
ETABSversion='v1' #2016 #no V16 #
NumIndiv=12;MutationPorc=20;NGenerations=25;LengthChromosome=0;elite=2;Nstopelite=6;NumberCores=6;tolbest=0.01;tolmean=0.01
#NumIndiv=16;MutationPorc=20;NGenerations=50;LengthChromosome=0;elite=2;Nstopelite=10;NumberCores=4

#import os;shareadress=os.getcwd()
Daughtercomputers=[]
#Daughtercomputers=[['asus',NumberCores,shareadress,NumIndiv]] #optional


##############################################################################################Process
import MAIN as MM
Pob_floatelite,fxelite,Pob_floatMean,Pob_floatmax=MM.etabsGA(Ke1Init,Element,ETABSversion,units,StructuralWeightObj,StructureDemand_CapacityObj,
                                                             ModelDirectory,ProgramPath,ModelName0,PerObj,ModalMassObj,NumIndiv,MutationPorc,
                                                             LengthChromosome,NGenerations,Nstopelite,elite,NumberCores,Daughtercomputers,tolbest,tolmean)