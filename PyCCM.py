# ~!/usr/bin/python
######################################################################
#date:2/27/2015   
#the pyccm 0.0 test program
#Auther: Junjie Jiang 
#Location: ASU Tempe AZ
#
#
#
#
######################################################################
#import module
import numpy as np
import scipy as sp
from scipy import spatial
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from math import *
#######################################################################     
def CCM(X,Y,X1,X2,X3):
	
	Xn=np.zeros((X2,X3))
	Yn=np.zeros((X2,X3))

	seed=np.random.rand(1,1)*(N-X3-X2)
	for i in range(0,X2):
		Xn[i,:]=X[int(seed)+i*X1:int(seed)+i*X1+X3]
		Yn[i,:]=Y[int(seed)+i*X1:int(seed)+i*X1+X3]
	Xn=np.transpose(Xn)
	Yn=np.transpose(Yn) 
	dist1=spatial.distance.cdist(Xn,Xn,'euclidean')
	dist2=spatial.distance.cdist(Yn,Yn,'euclidean')
	
	
	print (dist1.shape)

	#plt.imshow(dist1, cmap=plt.cm.jet)
	#plt.colorbar()
	#plt.show()

	#plt.imshow(dist2, cmap=plt.cm.jet) 
	#plt.colorbar()
	#plt.show()

	re_1=np.zeros((2,X3,X3))
	
#	print re_1.shape 
#	print dist2.shape
#	print Xn.shape
	re_1[0,:X3,:]=dist1
	re_1[1,:X3,:]=dist2	
#	re_1[0,(X3):(X3+X2), :]=np.transpose(Xn)
#	re_1[1,(X3):(X3+X2), :]=np.transpose(Yn)
#	print re_1
	number=np.zeros((X3,X2))
	for w in range(X3):
		#print w
		sit=re_1.shape
		temp3=np.zeros((sit[1],(X2+1)))
		temp3[:,0]=re_1[0,:,w]
		temp3[:,1:(X2+1)]=Yn
		#print '@@@'
		#print temp3
		temp3=temp3[temp3[:,0].argsort()]
		temp4=temp3[1:(X2+2),:]
		#print '&&&'
		#print temp3 
		sum1=0
		temp_1=np.zeros((X2+1,1))
		for po in range(X2+1):
			sum1=sum1+exp(-temp4[po,0]/temp4[0,0])
			temp_1[po,0]=exp(-temp4[po,0]/temp4[0,0])
		temp_1=temp_1/sum1
		#print '$$$'
		#print temp_1
		temp_1=np.transpose(temp_1)
		temp_2=temp4[:(X2+1),1:(X3+1)]
		#print temp_2
		#print '%%%'
		y_bar=np.dot(temp_1,temp_2)
		number[w,:]=y_bar;
	sum2=0
	for w1 in range(X2):
		sum2=sum2+pearsonr(number[:,w1],Yn[:,w1])[0]
	rho_y_x=sum2/X2
	
	print (rho_y_x)

	number=np.zeros((X3,X2))
	for w in range(X3):
		#print w
		sit=re_1.shape
		temp3=np.zeros((sit[1],(X2+1)))
		temp3[:,0]=re_1[1,:,w]
		temp3[:,1:(X2+1)]=Xn
		#print '@@@'
		#print temp3
		temp3=temp3[temp3[:,0].argsort()]
		temp4=temp3[1:(X2+2),:]
		#print '&&&'
		#print temp3 
		sum1=0
		temp_1=np.zeros((X2+1,1))
		for po in range(X2+1):
			sum1=sum1+exp(-temp4[po,0]/temp4[0,0])
			temp_1[po,0]=exp(-temp4[po,0]/temp4[0,0])
		temp_1=temp_1/sum1
		#print '$$$'
		#print temp_1
		temp_1=np.transpose(temp_1)
		temp_2=temp4[:(X2+1),1:(X3+1)]
		#print temp_2
		#print '%%%'
		x_bar=np.dot(temp_1,temp_2)
		number[w,:]=x_bar;
	sum2=0
	for w1 in range(X2):
		sum2=sum2+pearsonr(number[:,w1],Xn[:,w1])[0]
	rho_x_y=sum2/X2
	
	print (rho_x_y)

	return (rho_y_x,rho_x_y)	
######################################################################
#setting parameter
X1=int(1)
X2=int(2)

#X3=int(5000)
N=100000
#######################################################################
X=np.zeros(N)
Y=np.zeros(N)

X[0]=0.4
Y[0]=0.2
for i in range(1,N):
	X[i]=X[i-1]*(3.8-3.8*X[i-1]-0.02*Y[i-1])
	Y[i]=Y[i-1]*(3.5-3.5*Y[i-1]-0.1*X[i-1])
#########
#X is the first timeseries; Y is the second timeseries; X1 is the time lag;
#X2 is the embedded dimension; X3 is the used timeseries length
#########
rho_y_x=[]
rho_x_y=[]
lengt=[]
for i in range(10,3000,100):
	sum11=0
	sum22=0
	for j in range(10):
		rho_x, rho_y= CCM(X,Y,X1,X2,i)
		sum11=sum11+rho_x
		sum22=sum22+rho_y
	rho_y_x.append(sum11/10.0)
	rho_x_y.append(sum22/10.0)
	lengt.append(i)

plt.plot(lengt,rho_y_x,'bs--',lengt,rho_x_y,'g^--')
plt.show()



	
