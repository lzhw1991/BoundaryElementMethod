#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate


def prePar(fpar, distance, topography, moho, nrec):
    fp = open(fpar, 'w')
    context = \
        ['This is a parameter file for a 2-layer model with free-surface\n',
         '2 1 6 1 64.0 0.008 20 0.1 100\n',
         '2.0 1.0 0.0 3.0 1 1\n',
         '2 -60.0 840.0 3 5.0 90.0 0.0\n',
         str(nrec)+'\n',
         'init\n',
         'initu.dat\n',
         'initw.dat\n']
    nelement = len(distance)
    domain1 = ['domain 1 2.6 5.8 3.2 ' + str(2 * nelement) + '\n',
               '1 ' + str(nelement) + ' ' + str(nelement + 1)
               + ' ' + str(2 * nelement) + '\n']
    for line in context:
        fp.write(line)
    for line in domain1:
        fp.write(line)
    for i in range(nelement):
        fp.write('%d %f %f 0\n' % (i + 1, distance[-i-1], topography[-i-1]))
    for i in range(nelement):
        fp.write('%d %f %f 0\n' %
                 (i + 1 + nelement, distance[i], moho[i]))

    domain2 = ['domain 2 3.4 8.1 4.5 2\n',
               '1 2 3 4\n',
               '1 0.0 0.0 0\n',
               '2 0.0 0.0 0\n',
               '3 0.0 0.0 0\n',
               '4 0.0 0.0 0\n']
    for line in domain2:
        fp.write(line)
    fp.close()


def preReceivers(frec, distance, topography):
    fp = open(frec, 'w')
    nrec = len(distance)
    for i in range(nrec):
        fp.write('1 1 %f %f\n' % (distance[i], topography[i]))
    fp.close()


stations = np.load('mapped_stations.npy')

topography = stations[:, 2]
topography = np.append(topography[0], topography)
topography = np.append(topography, topography[-1])
topography = topography / 1e3

# interpTopography =
distance = stations[:, -1]
distance = np.append(distance[0] - 10e3, distance)
distance = np.append(distance, distance[-1] + 10e3)
distance = distance + 10e3
distance = distance / 1e3


# print distance

LAYER0 = 10  # top layer


plt.figure()
plt.plot(distance, topography, 'x')
plt.xlim([distance[0], distance[-1]])
# plt.ylim([0, 50])
# plt.axis('imag')

dxSrc = 1  # km
interpDistance = np.arange(distance[0], distance[-1], dxSrc)
#interpDistance = np.arange(200.0, 280, dxSrc)

f = interpolate.interp1d(distance, topography, kind='slinear')
interpTopography = f(interpDistance)
interpTopography = interpTopography

interpTopography = LAYER0 - interpTopography

print distance[0], distance[-1]

dist1 = np.arange(distance[0], 190)
dist2 = np.arange(220, distance[-1])
depth1 = dist1 * 0 + 40 + LAYER0
depth2 = dist2 * 0 + 70 + LAYER0
distMoho = np.append(dist1, dist2)
depthMoho = np.append(depth1, depth2)
fMoho = interpolate.interp1d(distMoho, depthMoho, kind='cubic')
interpMoho = fMoho(interpDistance)

prePar('ebem.dat', interpDistance, interpTopography, interpMoho, len(interpTopography))
preReceivers('receiver.dat', interpDistance, interpTopography)

plt.figure()
plt.plot(interpDistance, interpTopography, 'x')
plt.hold('True')
plt.plot(interpDistance, interpMoho)
plt.axis('image')
plt.xlim([interpDistance[0], interpDistance[-1]])
plt.ylim([LAYER0 - 10, 140])
plt.gca().invert_yaxis()

plt.show()
