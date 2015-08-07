#encoding: utf-8
from __future__ import division
import os
import sys
import string
import re
from numpy import*
from pulp import *
#from decimal import *


def dcc_lin(E,dim,mu):
    #getcontext().prec = 4
    ydim=mu[E[0,0]][E[0,1]].shape[2]
    n=dim.shape[0]

    #Find Neighbours of Node I
    N=[[] for i in range(n)]
    for e in range(E.shape[0]):
        i=E[e,0]
        j=E[e,1]
        N[i].append(j)
        N[j].append(i)
    

    #indicator function
    D=zeros((ydim,ydim))
    for y in range(ydim):
        for z in range(ydim):
            D[y,z]=-2*(y!=z)

    alpha=zeros((n,n))
    eta=zeros(n)
    beta=0

    #n*n dimension array
    i=0
    j=0
    nu=[[j for j in range(n)] for i in range(n)]
    delta=[[j for j in range(n)] for i in range(n)]
    gamma=[[j for j in range(n)] for i in range(n)]
    shift=0

    #nu:
    for e in range(E.shape[0]):
        i=E[e,0]
        j=E[e,1]

        #increase order step by one
        nu[i][j]=zeros((dim[i],dim[j],ydim))
        nu[i][j]=arange(dim[i]*dim[j]*ydim).reshape(dim[i],dim[j],ydim)+shift
        shift=shift+dim[i]*dim[j]*ydim
    #:alpha
    for e in range(E.shape[0]):
        i=E[e,0]
        j=E[e,1]
        alpha[i][j]=shift
        shift=shift+1
    #print 'alpha',alpha
    #eta:
    eta=arange(n)+shift
    #print 'eta',eta
    shift=shift+n
    #beta:
    beta=shift
    #print 'beta',beta
    shift=shift+1

    #delta:
    for e in range(E.shape[0]):
        i=E[e,0]
        j=E[e,1]
        delta[i][j]=zeros(dim[j])
        delta[i][j]=arange(dim[j])+shift
        shift=shift+dim[j]

        delta[j][i]=zeros(dim[i])
        delta[j][i]=arange(dim[i])+shift
        shift=shift+dim[i]
    
    #gamma:
    for e in range(E.shape[0]):
        i=E[e,0]
        j=E[e,1]
        gamma[i][j]=zeros((ydim,ydim))
        gamma[i][j]=arange(ydim*ydim).reshape(ydim,ydim)+shift
        shift=shift+ydim*ydim
    #print 'gamma',gamma

    num_vars=shift
    print "number of vars",shift

    #coeffience
    coef=zeros(num_vars)
    for e in range(E.shape[0]):
        i=E[e,0]
        j=E[e,1]
        for xi in range(dim[i]):
            for xj in range(dim[j]):
                for y in range(ydim):
                    coef[nu[i][j][xi,xj,y]]=mu[i][j][xi,xj,y]
    #print 'coef',coef
    

    #begin with LP

    #the LP variable:gives by a sequent
    
    VarSeq=range(num_vars)
    #print '@@@@@@@@@@@@@@@@@@@@@',VarSeq

    #set the problem as a min problem
    prob = LpProblem("TheDCCProblem", LpMinimize)

    #set the LP variable set a dictionary
    Var = LpVariable.dicts("VAR",VarSeq)
    #print 'Var',Var

    #set objection function
    prob += lpSum([coef[i]*Var[i] for i in VarSeq])
    #print 'object function',lpSum([coef[i]*Var[i]*1.0 for i in VarSeq])
    #set constraints
    #pulp.lpSum(vector):Calculate the sum of a list of linear expressions
    #Parameters: vector ¨C A list of linear expressions

    exprs=[]
    for e in range(E.shape[0]):
        i=E[e,0]
        j=E[e,1]
        exprs.append(Var[alpha[i][j]])
        
    i=0
    for i in range(n):
        exprs.append(Var[eta[i]])
        
    exprs.append(Var[beta])
    
    prob += lpSum(exprs) >= 0.0
    #print "lpSum",lpSum(exprs)

    
    for e in range(E.shape[0]):
        i=E[e,0]
        j=E[e,1]
        cur_alpha=alpha[i][j]
        for xi in range(dim[i]):
            for xj in range(dim[j]):
                for y in range(ydim):
                    for z in range(ydim):

                        if(y!=z):
                            prob += 1.0*Var[nu[i][j][xi,xj,y]]+1.0*Var[nu[i][j][xi,xj,z]]-1.0*Var[delta[i][j][xj]]-1.0*Var[delta[j][i][xi]]-1.0*Var[gamma[i][j][y,z]]-1.0*Var[cur_alpha] >= 0.0
                            #print 'constraint:',Var[nu[i][j][xi,xj,y]]+Var[nu[i][j][xi,xj,z]]-Var[delta[i][j][xj]]-Var[delta[j][i][xi]]-Var[gamma[i][j][y,z]]-Var[cur_alpha]


                        else:
                            prob += 1.0*Var[nu[i][j][xi,xj,y]]-1.0*Var[delta[i][j][xj]]-1.0*Var[delta[j][i][xi]]-1.0*Var[gamma[i][j][y,z]]-1.0*Var[cur_alpha] >= 0.0
                            #print 'constraint:',Var[nu[i][j][xi,xj,y]]-Var[delta[i][j][xj]]-Var[delta[j][i][xi]]-Var[gamma[i][j][y,z]]-Var[cur_alpha]
                        
    exprs=[]
    for i in range(n):
        cur_eta=eta[i]
        for xi in range(dim[i]):
            exprs=[]
           
            for j in N[i]:
                exprs.append(1.0*Var[delta[j][i][xi]])
            exprs.append(-1.0*Var[cur_eta])
            prob += lpSum(exprs)>=0.0
            #print 'delta constraint',lpSum(exprs)

    exprs=[]
    for y in range(ydim):
        for z in range(ydim):
            exprs=[]
            
            for e in range(E.shape[0]):
                i=E[e,0]
                j=E[e,1]
                exprs.append(1.0*Var[gamma[i][j][y,z]])
            exprs.append(-1.0*Var[beta])
            exprs.append(1.0*D[y,z])
           
            prob += lpSum(exprs)>=0.0
            #print 'gamma constraint',lpSum(exprs)


    prob.writeLP("The DCC Problem.lp")
    prob.solve()

    #save the result value
    res=[[j for j in range(n)] for i in range(n)]
    res1=[[j for j in range(n)] for i in range(n)]
    for e in range(E.shape[0]):
        i=E[e,0]
        j=E[e,1]
        res[i][j]=zeros((dim[i],dim[j],ydim))
        for xi in range(dim[i]):
            for xj in range(dim[j]):
                for y in range(ydim):
                    index=int(nu[i][j][xi,xj,y])
                    res[i][j][xi,xj,y]=Var[index].value()
                    #print '#',index,'   ','res:',res[i][j][xi,xj,y]
    #for v in prob.variables():
        #print v.name, "=", v.varValue
       
    
    return res

def PriorStat(E,dim,ydim,trainingPath):
    fileR = open(trainingPath)

    freq=[[j for j in range(n)] for i in range(n)]
    #init:
    for e in range(E.shape[0]):
        i=E[e,0]
        j=E[e,1]
        mu[i][j]=zeros((dim[i],dim[j],ydim))
        freq[i][j]=zeros((dim[i],dim[j],ydim))
    
    for line in fileR:
        line = line.split('\n')[0]
        line = line.split(',')
        for e in range(E.shape[0]):
            i=E[e,0]
            j=E[e,1]
            xi=line[i]
            xj=line[j]
            y=line[-1]
            freq[i][j][xi,xj,y]=freq[i][j][xi,xj,y]+1

    sum_count=0
    for e in range(E.shape[0]):
        i=E[e,0]
        j=E[e,1]
        sum_count=0

        for xi in range(dim[i]):
            for xj in range(dim[j]):
                for y in range(ydim):
                    sum_count=sum_count+freq[i][j][xi,xj,y]

        for xi in range(dim[i]):
            for xj in range(dim[j]):
                for y in range(ydim):
                    if(sum_count!=0):
                        mu[i][j][xi,xj,y]=freq[i][j][xi,xj,y]/sum_count                        

        

    fileR.close()    
    return mu


def TestModel(E,nu,testPath):
    W=0
    R=0
    print 'TEST MODEL##################'
    fileR = open(testPath)
    for line in fileR:
        line = line.split('\n')[0]
        line = line.split(',')
        sum_1=0
        sum_2=0
        for e in range(E.shape[0]):
            i=E[e,0]
            j=E[e,1]
            xi=line[i]
            xj=line[j]
            sum_1=sum_1+nu[i][j][xi,xj,0];
            sum_2=sum_2+nu[i][j][xi,xj,1];
            
           
        label=line[-1]
        #print '***********sum',sum_1,'   ',sum_2
        f1=-sum_1/2
        f2=-sum_2/2
        #print f1
        #print f2

        if(f1<f2):
            testlabel=1
        else:
            testlabel=0

        #print 'label',label
        #print 'testlabel',testlabel
        if(testlabel == int(label)):
            R=R+1
        else:
            W=W+1
    
    print 'R',R
    print 'W',W
    
            
if __name__ == "__main__":
   
    trainingPath="adulttraining.data"
    testPath="adulttest.test"
    #%E=[1,7;2,6;2,11;2,12;2,13;4,5;4,7;6,8;7,2;9,4;9,7;9,14;10,3;10,7;14,3;14,6];
    E=array([[0,6],[1,5],[1,10],[1,11],[1,12],[3,4],[3,6],[5,7],[6,1],[8,3],[8,6],[8,13],[9,2],[9,6],[13,2],[13,5]])
    dim=array([7,8,5,16,4,7,14,6,5,2,4,3,4,41])
    labelNum=2
    n=dim.shape[0]
    mu=[[j for j in range(n)] for i in range(n)]
    nu=[[j for j in range(n)] for i in range(n)]
    mu=PriorStat(E,dim,labelNum,trainingPath)
    nu=dcc_lin(E,dim,mu)
    TestModel(E,nu,testPath)
    print "done!"
    

          
