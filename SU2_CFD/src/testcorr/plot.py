import numpy as np, pylab as py

d = np.loadtxt('o')

py.figure()
py.subplot(2,2,1)
py.semilogy(d[:,0])
py.title('res')
py.subplot(2,2,2)
py.plot(d[:,1])
py.title('re_th')
py.subplot(2,2,3)
py.plot(d[:,2])
py.title('f_lamb')
py.subplot(2,2,4)
py.plot(d[:,3])
py.title('lamb')
py.show()
