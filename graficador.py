import pylab
import numpy
import math


data = numpy.loadtxt("trayectoria.dat")
t = data[:,0] 
Vel = data[:,1] 
Pos= data[:,2]
Costheta = data [:,3]
Sentheta= data [:, 4]
  
ejex= Pos*Costheta
ejey= Pos*Sentheta

pylab.plot(ejex,ejey, 'm')
pylab.xlabel('t')
pylab.ylabel('x')
pylab.title('X-T')
pylab.savefig("trayectoria.png")
pylab.close()
	
	
	
print "La grafica fue generada y guardada"
