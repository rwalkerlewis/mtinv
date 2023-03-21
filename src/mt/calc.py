import numpy.linalg
L = [ [1., -1., 0.], [-1., 2., -1.], [0., -1., 1.] ]
[eigenvalues,eigenvectors] = numpy.linalg.eig(L)
print L
print eigenvalues
print eigenvectors
quit()
