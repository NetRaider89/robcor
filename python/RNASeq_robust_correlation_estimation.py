import numpy as np
import scipy.io as sio
import scipy.sparse as sparse
import itertools as it
from cormad import cormad
from multiprocessing import JoinableQueue,Process

RNASeq = np.loadtxt("RNASeq.txt")
#RNASeq = np.loadtxt('rna_seq_full.txt', skiprows=1, usecols=range(1,152))

def Producer(inQueue):
    print('Producer start!')
    nGenes = RNASeq.shape[0]
    #fill upper triangular matrix of correlations
    #for i,j in it.product(range(nGenes),repeat=2):
    message = []
        
    for i in range(nGenes):
        if i % 100 == 0:
            print(i)
        for j in range(i+1,nGenes):
            if len(message)<100:
                message.append((i,j))
                #submit indices of genes i,j to queue        
            else:
                inQueue.put(message)
                message = []
    print('Producer put all tuple in the queue')
    print('Producer Done!')

def Worker(inQueue,outQueue):
    print('Worker start!')
    while True:
        message = inQueue.get()
        if message is None: 
            inQueue.task_done()
            print('Worker done!')
            return
        out_message = []
        for (i,j) in message:
            cor_val = cormad(RNASeq[i, :],RNASeq[j, :])
            out_message.append((i, j, cor_val))
        outQueue.put(out_message)
        inQueue.task_done()

def Consumer(outQueue):
    print('Consumer start!')
    CorMat = np.zeros((RNASeq.shape[0],RNASeq.shape[0]))
    while True:
        message = outQueue.get()
        if message is None:
            outQueue.task_done()
            print('Consumer done!')
            cor_mat_sp = sparse.csr_matrix(CorMat)
            sio.mmwrite('small_correlation_matrix.mtx', cor_mat_sp)
            return 
        for (i, j, cor_val) in message:
            CorMat[i, j]=cor_val
        outQueue.task_done()

def ParallelRobCor(X, n_jobs=24)
    inQueue = JoinableQueue(n_jobs)
    outQueue = JoinableQueue(n_jobs)
    
    producer_process = Process(target=Producer, args=(inQueue, n_jobs))
    worker_list = []
    for i in range(n_jobs):
       worker_list.append(Process(target=Worker, args=(inQueue, outQueue)))
    consumer_process =Process(target=Consumer, args=(outQueue,))
    producer_process.start()
    for w in worker_list:
        w.start()
    consumer_process.start()

    producer_process.join()
    for i in range(n_jobs):
        #Send poison message to processes
        inQueue.put(None)
    #quando arrivo a questo punto vuol dire che tutti i task 
    #nella coda di input sono finiti e il producer e i worker sono morti
    for w in worker_list:
        w.join()

    outQueue.put(None)
    consumer_process.join()
    #quando arrivo a questo punto tutti i task sulla coda di 
    #output sono stati smaltiti e anche il consumer muore

