# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 20:35:09 2022

@author: steve
"""

import numpy as np

def size(x):
    '''m=row n=column'''
    try:
        n=len(x[0][:])
        m=len(x)
    except TypeError:
        try:
            n=len(x)
            m=1
        except:
            n=1
            m=1
    return m,n
def rnum(i,j,k):
    num=k
    mat=[[num for a in range(1,j+1)]for a in range(1,i+1)]
    return mat    
def zeros(i,j):
    mat=rnum(i,j,0)
    return mat
def trans(x):
    try:
        c=max([len(x[i]) for i in range(len(x))])
        xt=zeros(c,len(x))
        for i in range(0,len(x)):
            for j in range(0,c):
                if len(x[i])>j:
                    xt[j][i]=x[i][j]
                else:xt[j][i]=0
        if c==1:
            xt=xt[0]
    except TypeError:
        try:
            xt=zeros(len(x),1)
            for i in range(0,len(x)):
                xt[i][0]=x[i]
        except TypeError:
            xt=x
    return xt

########################################################                  GA                      #####################################################################

def fPobDec(Pob_2,ne,Nr):
    Pob_10=[0]*Nr*ne
    Pob_2N=[]
    for i in range(Nr*ne):
        Pob_2Ni='0b'
        for j in range(len(Pob_2[0])):
            Pob_2Ni=Pob_2Ni+str(Pob_2[i][j])
        Pob_2N.append(Pob_2Ni)
        
    for i in range(Nr*ne):
        Pob_10[i]=int(Pob_2N[i],2)
    return Pob_10

def fPobDecInv(Pob_10,ne,Nr,LC):
    Pob_2=zeros(Nr*ne,LC)
    Pob_2N=[bin(Pob_10[i]) for i in range(len(Pob_10))]
    for i in range(Nr*ne):
        bi=Pob_2N[i]
        cont=LC
        for j in range(len(bi)-1,2-1,-1):
            cont-=1
            Pob_2[i][cont]=int(bi[j])
    return Pob_2

### NUMEROS BINARIOS - RANGO

def fScale(Pob_10,lim,ne,LC,Nr):
    Pob_float=[0]*Nr*ne
    for i in range(Nr):
        for j in range(ne):
            r=(i)*ne+j
            if lim[j][2]=='int':
                Pob_float[r]=int(round(Pob_10[r]/(2**LC-1)*(lim[j][1]-lim[j][0])+lim[j][0],0))
            elif lim[j][2]=='float':
                Pob_float[r]=float(round(Pob_10[r]/(2**LC-1)*(lim[j][1]*10000-lim[j][0]*10000)+lim[j][0]*10000,0)/10000)
    return Pob_float

def fScalInv(Pob_float,lim,ne,LC,Nr):
    Pob_10=[0]*Nr*ne
    for i in range(Nr):
        for j in range(ne):
            r=(i)*ne+j
            Pob_10[r]=int(round((Pob_float[r]-lim[j][0])*(2**LC-1)/(lim[j][1]-lim[j][0]),0))
    return Pob_10

###### PORCENTAJE EMPARENTARSE - RULETA - MAYOR AREA - PARAMETRO ALEATORIO  ######
def fSelectionRul(fxp,Nr):
    import math
    import numpy as np
    mean=sum(fxp)/len(fxp)
    if mean!=0:
        Ndecr=[fxp[i]/mean for i in range(len(fxp))]
        decr=[];resr0=[]
        for i in range(len(fxp)):
            decr.append(math.floor(Ndecr[i]))
            resr0.append(Ndecr[i]-decr[i])
        resr=[resr0[i]*360/sum(resr0) for i in range(len(resr0))]
        for i in range(1,len(resr)):
            resr[i]=resr[i]+resr[i-1] 
        while sum(decr)<Nr:
            roulette=360*np.random.rand()
            i=1
            while i<Nr-1:
                if roulette>=resr[i-1] and roulette<resr[i] and decr[i-1]<Nr/2:
                    decr[i-1]=decr[i-1]+1
                    i=Nr
                i+=1
    else:
        decr=[1]*len(fxp)
    return decr

###### FUNCION REPRODUCCION
###### MEZCLA DE GENES-CROMOSOMAS 

def freproduction(Pob_2,Nr,ne,decr,LC):
    import copy    
    import numpy as np
    dec=copy.deepcopy(decr)
    Parents=[0]*Nr
    cont=0
    while cont<Nr:
        v=max(dec)
        ind=dec.index(v)
        Parents[cont:cont+v]=[ind]*v
        dec[ind]=0
        cont+=v
    t1=np.random.permutation(Nr)
    if len(t1)%2!=0:print('The Number of Individuals must be even')
    t2=[[t1[i],t1[i+1]] for i in range(0,len(t1),2)]
    Couples=zeros(len(t2),2)
    for i in range(len(t2)):
        Couples[i][0]=Parents[t2[i][0]]
        Couples[i][1]=Parents[t2[i][1]]
    
    Pob_C=[0]*len(Pob_2)
    a1=[0]*len(t2)
    a2=[0]*len(t2)
    for t in range(len(t2)):
        a1[t]=Couples[t][0]
        a2[t]=Couples[t][1]
        for i in range(ne):
            pcross=np.random.randint(0,LC)
            r=(a1[t])*ne+i
            aux1=Pob_2[r][0:pcross]
            aux2=Pob_2[r][pcross:LC]
            r1=(a2[t])*ne+i
            aux3=Pob_2[r1][0:pcross]
            aux4=Pob_2[r1][pcross:LC]
            
            Pob_C[t*ne+i]=aux1+aux4
            Pob_C[int(len(Pob_2)/2)+t*ne+i]=aux3+aux2
    return Pob_C,Parents,Couples

##### MUTACION - CAMBIO AL AZAR

def fMutation(Pob_2,LC,MutationPorc):
    import numpy as np
    R,C=size(Pob_2)
    Tgenes=R*C
    GenesMut=round(Tgenes*MutationPorc/100)
    for i in range(GenesMut):
        row=np.random.randint(0,R)
        column=np.random.randint(0,C)
        Pob_2[row][column]=1-Pob_2[row][column]
    return Pob_2

##### FUNCION ELITE
##### QUEDAN LOS MEJORES ELEMENTOS DE TODAS LAS ITERACIONES
def fElite(fxp,fxpelite,elite,ne,Pob_2,punish,Pob_float,fx,Pob_2elite,fxelite,Pob_floatelite,multiobjectives):
    import copy
    import numpy as np
    if max(fxp)>fxpelite[0]:
        fxpd=copy.deepcopy(fxp)
        Pob_2elite0=[0]*(elite*ne)
        fxpelite0=[0]*(elite)
        for i in range(elite):
            indelite=fxpd.index(max(fxpd))
            fxpd[indelite]=0#statistics.mean(fxpd)
            ri=(indelite)*ne
            rj=(indelite+1)*ne
            Pob_2elite0[i*ne:(i+1)*ne]=Pob_2[ri:rj]
            fxpelite0[i]=fxp[indelite]
            if i==0 and punish[indelite]==0:
                Pob_floatelite=Pob_float[ri:rj]
                fxelite=fx[indelite]
                global contGelite
                contGelite=contG
        if multiobjectives==False or multiobjectives==[] or multiobjectives=='' or multiobjectives==0 or len(Pob_2elite)+1>=len(Pob_2)*0.3:
            Pob_2elite=copy.deepcopy(Pob_2elite0);fxpelite=copy.deepcopy(fxpelite0)#
        elif contG%2==0:Pob_2elite=copy.deepcopy(Pob_2elite0);fxpelite=copy.deepcopy(fxpelite0)
    for i in range(int(len(Pob_2elite)/ne)):
        indelite=int(round(np.random.rand()*(len(fxp)-1),0))
        fxp[indelite]=fxpelite[i]
        ri=(indelite)*ne
        rj=(indelite+1)*ne
        Pob_2[ri:rj]=Pob_2elite[i*ne:(i+1)*ne]
    return Pob_2,Pob_2elite,fxelite,Pob_floatelite,fxpelite,fxp


##### ALGORITMO GENETICO QUE UNE TODAS LAS DEMAS FUNCIONES
def GA(functionanalysis,otherdata,lim,NumIndiv,MutationPorc,LengthChromosome,NGenerations,Nstopelite,elite,\
       NumberCores,initialdata,Daughtercomputers,tolbest,tolmean,multiobjectives):
    
    import sys
    import os
    #### IMPORTACION LIBRERIA NUMPY
    if not 'numpy'in sys.modules.keys():
        import subprocess
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'numpy'])
    import copy
    import statistics
    import time
    
    try:
        ##### LIBRERIA QUE PERMITE TRABAJAR CON VARIOS NUCLEOS
        from joblib import Parallel, delayed
    except:
        if not 'joblib' in sys.modules.keys():
            try:
                import subprocess
                subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'joblib'])
                from joblib import Parallel, delayed
            except:
                exitl=0
                while exitl<=12:
                    try:
                        exitl+=1
                        if exitl==1:os.system("start powershell python -m pip install joblib")
                        time.sleep(15)
                        from joblib import Parallel, delayed
                        exitl=1000
                    except:
                        exitl+=1
    import multiprocessing
#    from sklearn.neighbors import KNeighborsClassifier
##### MACHINE LEARNING - NO SE UTILIZA
    from sklearn import linear_model
    clffx = linear_model.BayesianRidge()
    clfpunish = linear_model.BayesianRidge()
    ML=0
    import importlib
    global contG,stop
    
    tstart = time.time()
    Nr=NumIndiv;MutationPorc0=MutationPorc
    
       
    #### LIMITES - SE NECESITA ESTABLECER UN LIMITE - PAPER INCLUYE SUS LIMITES EN LA PROGRAMACION
    
    lim=[copy.deepcopy(lim[i]) for i in range(len(lim))]
    LCmin=max(5,len(bin(int(round(max([lim[i][1]-lim[i][0] for i in range(len(lim))]),0)+1)))-1)
    if LengthChromosome==[] or LengthChromosome=='':LengthChromosome=0
    if LCmin>LengthChromosome:print('LengthChromosome less than LengthChromosome required. LengthChromosome='+str(LCmin)+' is taken.')
    LC=max(LengthChromosome,LCmin)
    ne=len(lim)
    
#### NUMERO DE NUCLEOS CON LOS QUE TRABAJA - 4 ES OPTIMO
    num_cores = min(NumberCores,multiprocessing.cpu_count())
    
    
    functionanalysis.index('.');analysismodule='';aux=0;function=''
    for i in range(len(str(functionanalysis))):
        if functionanalysis[i]=='.':aux=1
        if aux==0:
            analysismodule=analysismodule+functionanalysis[i]
        elif functionanalysis[i]!='.':
           function=function+functionanalysis[i]
    if not os.path.exists("tempMultiComputer"):
        os.makedirs('tempMultiComputer')
    resultstotal='tempMultiComputer\\resultstotal'+str(int(round(np.random.rand()*100,0)))+time.strftime('_%d-%m-%y_%H-%M-%S')+'.out'
    with open(resultstotal,'w') as f0:
        f0.write('')
    
    
    #### SACAR RESULTADOS DE CODIGOS DE ETABS
    #### CREAR CARPETA TEMPORALY ALMACENAR INFO OBTENIDA DE OTROS PROGRAMAS
    
    resrand='Main'#str(int(round(np.random.rand()*100000,0)))
    if function!='py':
        Module = importlib.import_module(analysismodule)
        method=getattr(Module,function)
    else:
        from shutil import copyfile
        copyfile(functionanalysis, 'tempMultiComputer\\functionanalysiscopy1.py')
        with open('tempMultiComputer\\functionanalysiscopy1.py', 'a') as f:
            f.write('\r' + 'with open("tempMultiComputer\\'+'\\results'+resrand+'"+str(j0987654321)+".out", "a") as f:')
            f.write('\r' + '    res=str([x,fxi,punishi,individual0987654321]);res1=""')
            f.write('\r' + '    for i in range(len(res)):')
            f.write('\r' + '        if res[i]!=" ":res1=res1+res[i]')
            f.write('\r' + '    f.write("'+'\\'+'n"+res1)')
        data_list=open('tempMultiComputer\\functionanalysiscopy1.py', 'r').readlines()

        filessource=[arch.name for arch in os.scandir(os.getcwd()) if arch.is_file()]
        for j in range(len(filessource)):
            copyfile(filessource[j],'tempMultiComputer'+'\\'+filessource[j])
        
    num_coresall=[num_cores]
    createfun=[Nr]
    
    if function=='py':
        for j in range(0,num_coresall[0]):
            with open('tempMultiComputer\\functionanalysiscopy'+str(j)+'.py','w') as f:f.write('')
            
### INICIALIZA LA POBLACION - MATRICES DE CEROS

############Preparation of process
    Pob_2=zeros(Nr*ne,LC)

    Pob_floati=[]   #PARA CAMBIAR DE NUMEROS BINARIOS A ENTEROS - CREAR RANGOS
    if initialdata!=[]:
        for i in range(len(initialdata)):
            Pob_floati=Pob_floati+initialdata[i]
        Pob_10i=fScalInv(Pob_floati,lim,ne,LC,len(initialdata))
        Pob_2i=fPobDecInv(Pob_10i,ne,len(initialdata),LC)
        Pob_2[0:len(Pob_floati)]=Pob_2i
        
    for i in range(len(Pob_floati),Nr*ne):
        for j in range(LC):
                Pob_2[i][j]=np.random.randint(2)
    Pob_10= fPobDec(Pob_2,ne,Nr)
    Pob_float = fScale(Pob_10,lim,ne,LC,Nr)
    fxpelite=[0]*(elite)
    Pob_2elite=[];fxelite=100;Pob_floatelite=0;Pob_floatelite0=0
    contG=0;stop=0#;stopmean=0;Nstopmean=10#contelite=0    
    if tolbest==[] or tolbest=='':tolbest=0
    if tolmean==[] or tolmean=='':tolmean=0
    Mfxp=tolmean+1;meanfxp=0;meanfxp0=100;keepdata=[]
    if Nstopelite==0 or Nstopelite==[]: Nstopelite=NGenerations
    Pob_floatMean=0;Pob_floatMean0=100
    train_accuracy_clf=0
    
    
    #### COMIENZAN ALGORITMOS G
    #### CONDICIONES DE FINALIZACION
    while contG<NGenerations and stop<Nstopelite and fxelite>tolbest and abs(Mfxp-meanfxp)>tolmean \
        and abs(meanfxp-meanfxp0)>tolmean:# and Pob_floatMean0!=Pob_floatMean:
        contG+=1
        fxp=[]   ## FUNCIONES OBJETIVO  -- EN ESTE CASO  --PERIODOS Y MASA MODAL EFECTIVA
        fx=[]
        punish=[]
        if function=='py':
            for i in range(num_coresall[0]):
                ### CREA ARCHIVO DE RESULTADOS
                ### UTILIZA RESULTADOS PARA EVALUAR SI SE CUMPLEN LOS OBJETIVOS
                with open("tempMultiComputer\\results"+resrand+str(i)+".out", "w") as f2:
                    f2.write('[[0],0,0,"na"]')
                    f2.write('\r[[0],0,0,"na"]')
        
        
        ### ALANZAR PERIODO OBJETIVO Y MASA MODAL OBJETIVO - MULTIOBJETIVOS
        ### CREAR PONDERACION A OBJETIVOS 
        if isinstance(multiobjectives,int) and multiobjectives!=0:
            global weights
#            uk=[np.random.rand() for i in range(multiobjectives)];weights=[uk[i]/sum(uk) for i in range(len(uk))]
            if contG%2==0:  uk=[np.random.rand() for i in range(multiobjectives)];weights=[uk[i]/sum(uk) for i in range(len(uk))]
            else:           weights=[1/multiobjectives for i in range(multiobjectives)]
        else:weights=False
#        print(weights)
        

        #### FUNCION OBJETIVO
        #### DEPENDE DEL NUMERO DE NUCLEOS - PARA EVITAR SUPERPOSICION O ERROR DATOS
        #### XQ CORREN AL MISMO TIEMPO VARIOS PROGRAMAS
        def fMainAnalysis(i):
            import numpy as np
            import os
            import time
            j=i%num_cores
            ri=(i)*ne
            rj=(i+1)*ne
            x=Pob_float[ri:rj]
            if keepdata!=[] and trans(keepdata)[0].__contains__(x):
                indkeepdata=trans(keepdata)[0].index(x)
                fxi=keepdata[indkeepdata][1];punishi=keepdata[indkeepdata][2]
            else:
                if function=='py':
                    #### CARPETA TEMPORAL NOMBRA DE FORMA UNICA CADA PROGRAMA
                        data_list1=['x='+str(x)+'\n']+['individual0987654321='+str(i)+'\n']+['j0987654321='+str(j)+'\n']+data_list
                        fout=open('tempMultiComputer\\functionanalysiscopy'+str(j)+'.py', 'w')
                        fout.writelines(data_list1)
                        fout.close()
                        ### FUNCIONES Q AYUDAN POR SI EL PROGRAMA SE PARA O SUFRE ALGUN ERROR DURANTE
                        ### SU CORRIDA
                        os.system("start powershell python tempMultiComputer\\functionanalysiscopy"+str(j)+".py")
                        test=1
                        while test!=0:
                            try:
                                y=list(np.loadtxt('tempMultiComputer\\results'+resrand+str(j)+'.out',dtype='str'))
                                y=[eval(y[j]) for j in range(len(y))]
                                indy=trans(y)[3].index(i)
                                test=0
                            except:
                                time.sleep(5)
                        fxi=y[indy][1]
                        punishi=y[indy][2]
                else:
                    fxi,punishi=method(x,otherdata)
            return fxi,punishi ## RESULTADOS Q SE OBTIENEN - OBJETIVOS DE CASTIGO
        
        if contG>10 and ML!=0:
           for i in range(0,int(createfun[0])):
                ri=(i)*ne
                rj=(i+1)*ne
                x=Pob_float[ri:rj]
                fxi=clffx.predict(x)
                punishi=clfpunish.predict(x)
                results=[fxi,punishi]
                           ### MACHINE LAERNING - NO SE UTILIZA - ARRIBA
        else:
            if num_cores==1:
                results=[]
                for i in range(Nr):
                    resultsi=fMainAnalysis(i)
                    results.append(resultsi)
            else:
                if createfun[0]!=0:
                    #PERMITE TRABAJAR CON VARIOS NUCLEOS
                    results=Parallel(n_jobs=num_cores, backend="threading")(delayed(fMainAnalysis)(i) for i in range(0,int(createfun[0])))
    #                for i in range(0,int(createfun[0])):results=threading.Thread(target=fMainAnalysis, args=(i,));results.start()
                else:
                    results=[]

        if function=='py':
            with open('tempMultiComputer\\resultstemp.txt','w') as f: f.write(str(results))
            
        fx=trans(results)[0]
        punish=trans(results)[1]
        
        fxp=[]
        for i in range(len(fx)):
            K=punish[i]#*ne
            fxp.append(1/((fx[i]+K)*(K+1)))
        decr=fSelectionRul(fxp,Nr)
        
        ##################################################################Elite
        if contG>0: ## ALMACENA LOS IND ELITE
            Pob_2,Pob_2elite,fxelite,Pob_floatelite,fxpelite,fxp=\
            fElite(fxp,fxpelite,elite,ne,Pob_2,punish,Pob_float,fx,Pob_2elite,fxelite,Pob_floatelite,multiobjectives)
        Pob_floatMean0=Pob_floatMean; ##POBLACION PROMEDIO
        Pob_floatMean=[statistics.median(Pob_float[i:len(Pob_float):ne]) for i in range(ne)]
        printcont=1

        if function=='py':printcont=1

        ### DATOS Q SE GENERAN Y ALMACENAN EN UNA HOJA DE RESULTADOS
        meanfxp0=meanfxp
        meanfxp=statistics.mean(fxp);mfx=min(fx);Mfxp=max(fxp)
        indmax=fxp.index(Mfxp);ri=(indmax)*ne;rj=(indmax+1)*ne
        Pob_floatmax=Pob_float[ri:rj]
        if contG==1 or contG%printcont==0:
#            if Pob_floatelite!=0:
            sys.stdout.write('\r'+'elite Generation '+str(contGelite))
            sys.stdout.write(' elite fx '+str(round(fxelite,4))+ '  &  '+'elite Pob_float '+str(Pob_floatelite))
            sys.stdout.write('                                                                ')
            sys.stdout.write(' EndGeneration='+str(contG)+' minfx='+str(round(mfx,4))+' Pob_floatmax='+str(Pob_floatmax)+'               ')
            sys.stdout.flush()
         
        with open(resultstotal,'a') as f0:
            mpunish=min(punish);
            if contG==1:
                f0.write('\rGeneration0_meanfxp1_maxfxp2_minfx3_minpunish4_fxelite5_Pobfloatelite6_nameoffunctionanalized7=[')
            f0.write('\r['+str(contG)+','+str(meanfxp)+','+str(Mfxp)+','+str(mfx)+',')
            f0.write(str(mpunish)+','+str(fxelite)+','+str(Pob_floatelite))
            f0.write(',"'+functionanalysis+'"]')
        for jj in range(len(fx)):
            if keepdata!=[] and trans(keepdata)[0].__conatains__(Pob_float[(jj)*ne:(jj+1)*ne])==0 and ML!=0:
                keepdata.append([[[Pob_float[(jj)*ne:(jj+1)*ne]],fx[jj],punish[jj]]])
                if contG>5:
                    clffx.fit(trans(keepdata)[0],trans(keepdata)[1])
                    clfpunish.fit(trans(keepdata)[0],trans(keepdata)[2])
                    train_accuracy_clf = clffx.score(trans(keepdata)[0],trans(keepdata)[1])
                    if train_accuracy_clf>0.97:
                        ML+=1;print("train_accuracy_clf>",0.97)
                
        ##################################################################Selection
        Pob_2,Parents,Couples=freproduction( Pob_2,Nr,ne,decr,LC)
        Pob_2=fMutation(Pob_2,LC,MutationPorc)
        Pob_10= fPobDec(Pob_2,ne,Nr)
        Pob_float = fScale(Pob_10,lim,ne,LC,Nr)
        if contG==1:rewrite='w'
        else:rewrite='a'
        with open('tempMultiComputer\\dataGA.txt',rewrite) as f0:
            f0.write('\rfxp='+str(fxp))
            f0.write('\rdecr='+str(decr))
            f0.write('\rParents='+str(Parents))
            f0.write('\rCouples='+str(Couples))
            f0.write('\rPob_floatMean='+str(Pob_floatMean)+'\n')
        if Pob_floatelite!=0 and Pob_floatelite0==Pob_floatelite:
            stop+=1
        else:
            stop=0;MutationPorc=MutationPorc0
        if stop==round(Nstopelite*0.7,0):
            MutationPorc=MutationPorc0*2
        Pob_floatelite0=copy.deepcopy(Pob_floatelite)


    #### CONDICIONES PARA Q EL PROGRAMA PARE Y SUS AVISOS
    if stop>=Nstopelite:print('stop>=Nstopelite')
    if fxelite<=tolbest:print('fxelite<=tolbest')
    if abs(Mfxp-meanfxp)<=tolmean:print('abs(Mfxp-meanfxp)<=tolmean')
    if abs(meanfxp-meanfxp0)<=tolmean:print('abs(meanfxp-meanfxp0)<=tolmean')
    if Pob_floatMean0==Pob_floatMean:print('Pob_floatMean0==Pob_floatMean')
    
    
    with open(resultstotal,'a') as f0:
        f0.write('\r]')
    if function=='py':
        for j in range(0,num_coresall[0]):
            os.remove('tempMultiComputer\\functionanalysiscopy'+str(j)+'.py')
            os.remove('tempMultiComputer\\results'+resrand+str(j)+'.out')

    elapsed = time.time() - tstart
    print('\n elapsed [s]=',elapsed)
    drawresultsGA(resultstotal)
    return Pob_floatelite,fxelite,Pob_floatMean,Pob_floatmax


### FUNCION DIBUJAR RESULTADOS

def drawresultsGA(resultstotal):
    f=open(resultstotal,'r')
    resultstotal_list0=f.readlines()
    resultstotal_list=[]
    for i in range(2,len(resultstotal_list0)):
        try:resultstotal_list.append(eval(resultstotal_list0[i]))
        except:pass
    fxmeandraw=[];fxMaxdraw=[]
    for i in range(len(resultstotal_list)):
        fxmeandraw.append(1/resultstotal_list[i][1])
        fxMaxdraw.append(1/resultstotal_list[i][2])
    import matplotlib.pyplot as plt
    steps=list(range(len(resultstotal_list)))
    plt.plot(steps,fxmeandraw,c='silver',label='f(x) mean')
    plt.plot(steps,fxMaxdraw,c='black',label='f(x) best')
    plt.ylabel('f(x)');plt.xlabel('step')
    plt.legend()
    
    
########################################################              ETABS GA                   #####################################################################
###########################################################################################################################################################################

def etabselem(KeNew,Element,units,ModelDirectory,ProgramPath,ModelName0,ETABSversion,*remove):
    import numpy as np
    import comtypes.client
    import sys
    from shutil import copyfile
    import os
    
    Periods=[];Ux=[];Uy=[];Rz=[];contp=0
    while Periods==[] and Ux==[] and Uy==[] and Rz==[] and contp<3:
            contp+=1
            #try:
            #set the following flag to True to manually specify the path to ETABS.exe
            #this allows for a connection to a version of ETABS other than the latest installation
            #otherwise the latest installed version of ETABS will be launched
            if ProgramPath==[]:SpecifyPath=False
            else:SpecifyPath=True
            ### CREA NOMBRE ARCHIVOS - EVITAR MISMO NOMBRE PARA DOS ARCHIVOS
            ModelName=ModelName0+str(round(np.random.rand()*10e14))
            ModelPath0=ModelDirectory+ModelName0+'.EDB'
            ModelPath=ModelDirectory+ModelName+'.EDB'
            copyfile(ModelPath0, ModelPath)
            ### SE MANTIENE ARCHIVO ORIGINAL Y SE CREAN COPIAS TEMPORALES
            
            #set the following flag to True to attach to an existing instance of the program
            #otherwise a new instance of the program will be started
            AttachToInstance = False
            if AttachToInstance:
            #attach to a running instance of ETABS
                try:
                    #get the active ETABS object
                    myETABSObject = comtypes.client.GetActiveObject("CSI.ETABS.API.ETABSObject") 
                except (OSError, comtypes.COMError):
                    print("No running instance of the program found or failed to attach.")
                    sys.exit(-1)
            else:
                #create API helper object
                contdel=0
                while contdel<3:
                    try:
                        helper = eval("comtypes.client.CreateObject('ETABS"+ETABSversion+".Helper')")
                        helper = eval("helper.QueryInterface(comtypes.gen.ETABS"+ETABSversion+".cHelper)")
                        contdel=3
                        if SpecifyPath:
                            try:
                                #'create an instance of the ETABS object from the specified path
                                myETABSObject = helper.CreateObject(ProgramPath)
                            except (OSError, comtypes.COMError):
                                print("Cannot start a new instance of the program from " + ProgramPath)
                                sys.exit(-1)
                        else:
                    
                            try: 
                                #create an instance of the ETABS object from the latest installed ETABS
                                myETABSObject = helper.CreateObjectProgID("CSI.ETABS.API.ETABSObject") 
                            except (OSError, comtypes.COMError):
                                print("Cannot start a new instance of the program.")
                                sys.exit(-1)  
                    except:contdel+=1
            
            try:
                #start ETABS application
                myETABSObject.ApplicationStart()
            #create SapModel object
                SapModel=myETABSObject.SapModel
            
            #initialize model   
                SapModel.InitializeNewModel()
                #create file
                ### MODEL PATH - NOMBRE DEL ARCHIVO
                ### ASSERT - EVALUA RESULTADOS 
                assert SapModel.File.OpenFile(ModelPath)==0,'error '+ModelPath
                ret = SapModel.File.OpenFile(ModelPath)
#                print(ModelPath,ret)
                assert SapModel.SetModelIsLocked(False)==0,'error '+ModelPath
                ret = SapModel.SetModelIsLocked(False)
                
#               """switch to N-m units"""
                unit=14 #"""kgf_cm_C"""
                assert SapModel.SetPresentUnits(unit)==0,'error '+ModelPath
                SapModel.SetPresentUnits(unit)
                for j in range(len(Element)):
                    for k in range(len(Element[j])-1):
                        if Element[j][0]=='Link1':
                            Name=str(Element[j][k+1])
                            DOF=[1, 0, 0, 0, 0, 0]
                            Fixed=[0, 0, 0, 0, 0, 0]
                            NonLinear = [1,0,0,0,0,0]
                            Keq=[KeNew[(j+1)*k], 0, 0, 0, 0, 0]
                            Yield = [0.5633*KeNew[(j+1)*k],0,0,0,0,0]
                            Keff=[0.09091*KeNew[(j+1)*k], 0, 0, 0, 0, 0]
                            Ce=[0, 0, 0, 0, 0, 0] #AMORTIGUAMIENTO
                            Ratio = [0.023,0,0,0,0,0]
                            Exp = [1.2,0,0,0,0,0]
                            print('======1======')
                            #SapModel.PropLink.SetLinear(Name, DOF, Fixed, Ke, Ce,0,0) 
                            SapModel.PropLink.SetPlasticWen(Name, DOF, Fixed, NonLinear, Keq, Ce, Keff, Yield, Ratio, Exp, 0, 0)
                
                #SapModel.Analyze.SetRunCaseFlag("MODAL", False)
                #======= Linea guardar datos 
                #with open('rigidez.txt','a') as gg: gg.write('Ke='+str(KeNew)+'\r')

                assert SapModel.Analyze.RunAnalysis()==0,'error '+ModelPath
                SapModel.Analyze.RunAnalysis()

                ### DATOS A COMPARAR CON OBJETIVO - LUEGO SE BORRAN
                ### Extraccion derivas
                SapModel.Results.Setup.DeselectAllCasesAndCombosForOutput()
                assert SapModel.Results.Setup.SetCaseSelectedForOutput('953-1')==0,'error'+ModelPath
                SapModel.Results.Setup.SetCaseSelectedForOutput('953-1')

                #DERIVAS
                ResModal= SapModel.Results.StoryDrifts()
                Periods=ResModal[6] #Derivas
                
                # separar maximos en x ^ y
                # derivas=[x5 y5 x4 y4 ... x1 y1]
                der = np.zeros(10)
                #componentes x
                der[0]=Periods[0]
                der[2]=Periods[4]
                der[4]=Periods[8]
                der[6]=Periods[12]
                der[8]=Periods[16]
                #componentes y
                der[1]=Periods[1]
                der[3]=Periods[5]
                der[5]=Periods[9]
                der[7]=Periods[13]
                der[9]=Periods[17]
                
                Periods = der
                
                print('COMPLETAOOOOO')
                
                ResModalParticipation=SapModel.Results.ModalParticipatingMassRatios()
                Ux=ResModalParticipation[5]
                Uy=ResModalParticipation[6]
                Rz=ResModalParticipation[13]
                
                # Unlock model
                if remove!=():
                    rem=remove[0]
                else:rem='ok'
                
#                if min(rem!='no',rem!='No',rem!='NO'):
                if rem!='no' or rem!='No' or rem!='NO':
                    print('API CERRRAAARR')
                    assert SapModel.SetModelIsLocked(False)==0
                    SapModel.SetModelIsLocked(False)
                    assert myETABSObject.ApplicationExit(False)==0
                    myETABSObject.ApplicationExit(False)
                    os.remove(ModelPath)
                    os.remove(ModelDirectory+ModelName+'.ico')
                    os.remove(ModelDirectory+ModelName+'.$et')
                    os.remove(ModelDirectory+ModelName+'.ebk')
                    
            # SI HAY ERROR AL INICIO - BORRA TODO ETABS
            except:
                print('REMOVEEEER')
                os.remove(ModelPath)
                try:os.remove(ModelDirectory+ModelName+'.ico')
                except:pass
                try:os.remove(ModelDirectory+ModelName+'.$et')
                except:pass
                try:os.remove(ModelDirectory+ModelName+'.ebk')
                except:pass
    return Periods,Ux,Uy,Rz
    #return der,Ux,Uy,Rz

# =============================================================================
# 
# =============================================================================

#except:
#    pass
def etabsIterations(KeNew,otherdata):
    import time
    Element=otherdata[0];units=otherdata[1];StructuralWeigthObj=otherdata[2];StructureDemand_CapacityObj=otherdata[3]
    ModelDirectory=otherdata[4];ProgramPath=otherdata[5];ModelName0=otherdata[6]
    PerObj=otherdata[7];ModalMassObj=otherdata[8];ETABSversion=otherdata[9]
    ob=len(PerObj)
    
###""" REESTRUCTURAR PORQUE SOLO FUNCIONA PARA LOS PERIODOS Y MASA MODAL EFECTIVA"""

    Periods,Ux,Uy,Rz=etabselem(KeNew,Element,units,ModelDirectory,ProgramPath,ModelName0,ETABSversion)
    ## MASA MODAL EFECTIVA PARA LOS DIFRENTES MODOS DE VIBRACION
    while ob>len(Periods):
        try:
            time.sleep(3)
            Periods,Ux,Uy,Rz=etabselem(KeNew,Element,units,ModelDirectory,ProgramPath,ModelName0,ETABSversion)
        except:pass
        
        ### MODIFICA EL PERIODO PARA NO TOMAR EN CUENTA SUS RESULTADOS
    if ob>len(Periods):
        print('the number of objectives is greater than the number of results')
        error=2 ; punish=1 ; Periods=[100 for i in range(ob)]
    else:
        ## PONDERACIONES A OBJETIVOS
        vectOp=PerObj[:ob]
        vectOm=[]
        if ModalMassObj!=[]:
            for i in range(len(ModalMassObj)):vectOm=vectOm+trans(ModalMassObj)[i]
            Wp=0.8;Wm=0.2
        else:Wp=1;errorm=0
        vectp=Periods[:ob]
        vectm=Ux[:ob]+Uy[:ob]+Rz[:ob]
#        vect=vectp+vectm
        if weights!=False or weights!=[]:ukp=weights[:ob]
        else:ukp=list(range(len(vectp),0,-1))
        wkp=[ukp[i]/sum(ukp) for i in range(len(ukp))]
#        ukm=[np.random.rand() for i in range(len(vectm))]
        if weights!=False or weights!=[]:ukm=weights[ob:]
        else:ukm=list(range(len(vectm),0,-1))
        wkm=[ukm[i]/sum(ukm) for i in range(len(ukm))]
        
        errorp=sum([abs((vectp[i]-vectOp[i])/vectOp[i])*wkp[i] if vectOp[i]!=0 else vectp[i]*wkp[i] for i in range(len(vectp))])*Wp
        if ModalMassObj!=[]:errorm=sum([max(0,abs((vectm[i]-vectOm[i])/vectOm[i])-0.25)*wkm[i] if vectOm[i]!=0 else max(0,vectm[i]-0.25)*wkm[i] for i in range(len(vectm))])*Wm
        error=errorp+errorm
        punish=0
        with open('test.txt','w') as ff: ff.write('weights='+str(weights)+'\rwkp='+str(wkp)+' wkm='+str(wkm)+'\r vectOp='+str(vectOp)+' vectp='+str(vectp)+'\rvectOm='+str(vectOm)+' vectm='+str(vectm)+'\rerrorp='+str(errorp)+' errorm='+str(errorm))
    return error,punish


def etabsGA(Ke1Init,Element,ETABSversion,units,StructuralWeigthObj,StructureDemand_CapacityObj,ModelDirectory,
            ProgramPath,ModelName0,PerObj,ModalMassObj,NumIndiv,MutationPorc,LengthChromosome,NGenerations,
            Nstopelite,elite,NumberCores,Daughtercomputers,tolbest,tolmean):
    lim=[]
    initialdata=[]
    for i in range(len(Element)):
        if Element[i][0]=='Link1':
            for j in range(0,len(Element[i])-1):
                if Ke1Init!=[]:
                    # PROPIOS LIMITES - PROBAR NUMEROS PARA DETERMINAR Q FUNCIONA CORRECTAMENTE
                    #lim.append([Ke1Init[i][j]/9,Ke1Init[i][j]*9,'float'])
                    lim.append([Ke1Init[i][j],Ke1Init[i][j]*5,'float'])
                    initialdata.append(Ke1Init[i][j])
                else:
                    lim.append([0,1e7,'float'])
    otherdata=[Element,units,StructuralWeigthObj,StructureDemand_CapacityObj,ModelDirectory,ProgramPath,ModelName0,PerObj,ModalMassObj,ETABSversion]#optional
    functionanalysis='MAIN.etabsIterations'
    initialdata=[initialdata]
    if NumberCores>1:controlPowershell('ETABS',NumberCores)
        
    
    multiobjectives=len(PerObj)+sum([len(ModalMassObj[jj]) for jj in range(len(ModalMassObj))])
#    multiobjectives=0
    
        
    Pob_floatelite,fxelite,Pob_floatMean,Pob_floatmax=GA(functionanalysis,otherdata,lim,NumIndiv,\
            MutationPorc,LengthChromosome,NGenerations,Nstopelite,elite,NumberCores,initialdata,Daughtercomputers,tolbest,tolmean,multiobjectives)
    return Pob_floatelite,fxelite,Pob_floatMean,Pob_floatmax


###ERRORES ANORMALES DETENIAN LOS PROCESOS
### IDENTIFICA ANORMALIDADES PARA CORRER DONDE SE DETUVO
def controlPowershell(program,NumberCores):
    text=[]
    text.append("import os")
    text.append("\rimport psutil")
    text.append("\rimport time")
    text.append("\rimport copy")
    text.append("\rimport numpy as np")
    text.append("\rprint('Etabs Monitoring')")
    text.append("\rtol=5")
    text.append("\rcontdel=0;maxcontsaved=0")
    text.append("\rif '"+program+"'=='ETABS' or '"+program+"'=='Etabs' or '"+program+"'=='etabs':program1='ETABS.exe';program2='etabs.exe'")
    text.append("\rreppid=[];cont1=0;maxcont=0;pidetabs0=[];switch=0;cont=0")
    text.append("\rwhile True:")
    text.append("\r    pidetabs=[]")
    text.append("\r    for proc in psutil.process_iter(attrs=['pid', 'name']):")
    text.append("\r        if proc.info['name']==program1:")
    text.append("\r            pidetabs.append(proc.info['pid'])")
    text.append("\r            if reppid==[] or list(np.transpose(reppid)[0]).__contains__(proc.info['pid'])==0:")
    text.append("\r                reppid.append([proc.info['pid'],0]);cont1=0;switch=0;contdel=0")
    text.append("\r    if cont1==0:cont+=1;maxcontsaved=max(maxcontsaved,maxcont)")
    text.append("\r    cont1+=1")
    text.append("\r    if pidetabs0!=[] and pidetabs!=[]:")
    text.append("\r        for i in range(len(pidetabs)):")
    text.append("\r            for j in range(len(pidetabs0)):")
    text.append("\r                ind=list(np.transpose(reppid)[0]).index(pidetabs[i])")
    text.append("\r                if (pidetabs[i]==pidetabs0[j] and len(pidetabs)<=max(2,"+str(NumberCores)+"/1.5)) or (cont>"+str(NumberCores)+"*1.5 and cont1/3>maxcontsaved):")
    text.append("\r                    reppid[ind][1]=reppid[ind][1]+1.0")#+1")print(cont,maxcontsaved,reppid[ind]);
    text.append("\r                    if (min(list(np.transpose(reppid)[1]))>tol and cont1>maxcont) or cont1/3>maxcontsaved:")
    text.append("\r                        if switch==0:switch=1;cont1=0;reppid=[[reppid[j][0],0] for j in range(len(reppid))]")
    text.append("\r                        if (switch==1 and cont1>maxcont) or cont1/3>maxcontsaved:")
    text.append("\r                            os.system('taskkill /f /im ' + program2)")
    text.append("\r                            reppid[ind][1]=0;cont1=0;switch=0;cont=0;maxcont=maxcont/1.5;maxcontsaved=maxcontsaved/1.5")
    text.append("\r        for i in range(len(reppid)-1,-1,-1):")
    text.append("\r            if pidetabs.__contains__(reppid[i][0])==0:")
    text.append("\r                reppid.remove(reppid[i])")
    text.append("\r    pidetabs0=copy.deepcopy(pidetabs)")
    text.append("\r    print(reppid,cont1)")
    text.append("\r    time.sleep(5);maxcont=max(cont1,maxcont)")
    text.append("\r    if len(pidetabs)==0:")
    text.append("\r        contdel+=1")
    text.append("\r    if contdel==5:")
    text.append("\r       break")
    import os
    Daughtercomputers = os.getcwd()
    print(Daughtercomputers)
    with open(Daughtercomputers + '\controlGAprocess.py','w+') as f:f.writelines(text)
   # with open(Daughtercomputers + '\\tempMultiComputer\\controlGAprocess.py','w+') as f:f.writelines(text)
   # os.system("start powershell python tempMultiComputer\\controlGAprocess.py")
    os.system("start powershell python \controlGAprocess.py")
