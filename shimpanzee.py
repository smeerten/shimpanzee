#!/usr/bin/env python

# Copyright 2016 - 2017 Bas van Meerten and Wouter Franssen

# This file is part of Shimpanzee.
#
# Shimpanzee is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Shimpanzee is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Shimpanzee. If not, see <http://www.gnu.org/licenses/>.

import sip
import sys
import os
sip.setapi('QString', 2)
try:
    from PyQt4 import QtGui, QtCore
    from PyQt4 import QtGui as QtWidgets
    QT = 4
except ImportError:
    from PyQt5 import QtGui, QtCore, QtWidgets
    QT = 5
    
if __name__ == '__main__':
    root = QtWidgets.QApplication(sys.argv)
    root.setWindowIcon(QtGui.QIcon(os.path.dirname(os.path.realpath(__file__)) + '/logo.gif'))
    splash_pix = QtGui.QPixmap(os.path.dirname(os.path.realpath(__file__)) + '/logo.gif')
    splash = QtWidgets.QSplashScreen(splash_pix, QtCore.Qt.WindowStaysOnTopHint)
    splash.setMask(splash_pix.mask())
    progressBar = QtWidgets.QProgressBar(splash)
    progressBar.setGeometry(2.5*splash.width()/10, 0.96*splash.height(),5*splash.width()/10, 0.04 * splash.height())
    splash.show()    
    
splashSteps=8.0/100
splashStep = 0.0
def splashProgressStep(splashStep): #A function to easily increase the progressbar value
    if __name__ == '__main__':
        splashStep=splashStep+1
        progressBar.setValue(splashStep // splashSteps + (splashStep % splashSteps > 0)) #Rounds up without math or numpy module
        root.processEvents()   
    return splashStep


    
import matplotlib
splashStep = splashProgressStep(splashStep)
if QT ==4:
    matplotlib.use('Qt4Agg')
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
else:
    matplotlib.use('Qt5Agg')
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
splashStep = splashProgressStep(splashStep)
from matplotlib.figure import Figure
splashStep = splashProgressStep(splashStep)
import numpy as np
splashStep = splashProgressStep(splashStep)
import matplotlib.pyplot as plt
splashStep = splashProgressStep(splashStep)
import scipy.stats as stats
splashStep = splashProgressStep(splashStep)
from safeEval import safeEval
splashStep = splashProgressStep(splashStep)
from spectrumFrame_Shim import Plot1DFrame
splashStep = splashProgressStep(splashStep)
from scipy.interpolate import UnivariateSpline
splashStep = splashProgressStep(splashStep)



                    # x  y  z 
convert_1 = np.array([[0, 1.0, 0], 
                     [0, 0, 1], 
                     [1, 0, 0]])
rot_1 = np.array([[1, 0, 0], # x
                  [0, 0, -1], # y
                  [0, 1, 0]]) # z
convert_1_inv = np.dot(convert_1.T,np.linalg.inv(np.dot(convert_1,convert_1.T))) # Calculate the right inverse
X90_1 = np.dot(np.dot(convert_1, rot_1), convert_1_inv)
X_90_1 = np.linalg.inv(X90_1)

                     # xx xy xz yy yz zz
convert_2 = np.array([[0, 1, 0, 0, 0, 0],
                      [0, 0, 0, 0, 1, 0],
                      [-0.5, 0, 0, -0.5, 0, 1],
                      [0, 0, 1, 0, 0, 0],
                      [0.5, 0, 0, -0.5, 0, 0]])
convert_2[2] *= np.sqrt(1/3.0)
rot_2 = np.array([[1, 0, 0, 0, 0, 0], # xx
                  [0, 0, -1, 0, 0, 0], # xy
                  [0, 1, 0, 0, 0, 0], # xz
                  [0, 0, 0, 0, 0, 1], # yy
                  [0, 0, 0, 0, -1, 0], # yz
                  [0, 0, 0, 1, 0, 0]]) # zz
convert_2_inv = np.dot(convert_2.T,np.linalg.inv(np.dot(convert_2,convert_2.T))) # Calculate the right inverse
X90_2 = np.dot(np.dot(convert_2, rot_2), convert_2_inv)
X_90_2 = np.linalg.inv(X90_2)

                     # xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz
convert_3 = np.array([[0, 3.0, 0, 0, 0, 0, -1, 0, 0, 0],
                      [0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                      [0, -1, 0, 0, 0, 0, -1, 0, 4, 0],
                      [0, 0, -3, 0, 0, 0, 0, -3, 0, 2],
                      [-1, 0, 0, -1, 0, 4, 0, 0, 0, 0],
                      [0, 0, 1, 0, 0, 0, 0, -1, 0, 0],
                      [1, 0, 0, -3, 0, 0, 0, 0, 0, 0]])
convert_3[0] *= 0.25*np.sqrt(35.0/2.0)
convert_3[1] *= 0.5*np.sqrt(105)
convert_3[2] *= 0.25*np.sqrt(21.0/2.0)
convert_3[3] *= 0.25*np.sqrt(7)
convert_3[4] *= 0.25*np.sqrt(21.0/2.0)
convert_3[5] *= 0.25*np.sqrt(105)
convert_3[6] *= 0.25*np.sqrt(35.0/2.0)
rot_3 = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0], # xxx
                  [0, 0, -1, 0, 0, 0, 0, 0, 0, 0], # xxy
                  [0, 1, 0, 0, 0, 0, 0, 0, 0, 0], # xxz
                  [0, 0, 0, 0, 0, 1, 0, 0, 0, 0], # xyy
                  [0, 0, 0, 0, -1, 0, 0, 0, 0, 0], # xyz
                  [0, 0, 0, 1, 0, 0, 0, 0, 0, 0], # xzz
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, -1], # yyy
                  [0, 0, 0, 0, 0, 0, 0, 0, 1, 0], # yyz
                  [0, 0, 0, 0, 0, 0, 0, -1, 0, 0], # yzz
                  [0, 0, 0, 0, 0, 0, 1, 0, 0, 0]]) # zzz
convert_3_inv = np.dot(convert_3.T,np.linalg.inv(np.dot(convert_3,convert_3.T))) # Calculate the right inverse
X90_3 = np.dot(np.dot(convert_3, rot_3), convert_3_inv)
X_90_3 = np.linalg.inv(X90_3)

Z1LIM = 20.0
X1LIM = 20.0
Y1LIM = 20.0

Z2LIM = 5.0
XZLIM = 10.0
XYLIM = 10.0
YZLIM = 10.0
X2_Y2LIM = 10.0

Z3LIM = 1.0
XZ2LIM = 1.0
YZ2LIM = 1.0
ZX2_ZY2LIM = 1.0
XYZLIM = 1.0
X3LIM = 1.0
Y3LIM = 1.0

NSTEPS = 1000

def rotMatrix_z(size, angle):
    return np.diag(np.cos(angle*np.arange(size,-size-1,-1)))+np.fliplr(np.diag(np.sin(angle*np.arange(size,-size-1,-1))))


class ShimSim(object):
    
    def __init__(self, N=100000, dimensions=[2.5,10], angles=[0,0], sw=500.0, npoints=512):
        super(ShimSim, self).__init__()
        self.OVERSAMPLE = 16
        self.N = N
        self.r = dimensions[0] # Diameter to radius [mm]
        self.l = dimensions[1] # length [mm]
        self.theta = angles[0]
        self.phi = angles[1]
        self.sw = sw # sweepwidth [kHz]
        self.npoints = npoints # number of points on the frequency axis
        self.fidSurface = 0.0
        self.fwhm = 0.0
        self.resetGame()
        self.setupGrid()

    def setupGrid(self):
        self.z = np.random.uniform(-self.l/2.0, self.l/2.0, self.N)
        R = np.random.triangular(0, self.r, self.r, self.N)
        angle = np.random.uniform(0, 2*np.pi, self.N)
        self.x = R*np.cos(angle)
        self.y = R*np.sin(angle)
        self.rotMatrix_1 = np.dot(np.dot(X_90_1, rotMatrix_z(1, self.theta)), np.dot(X90_1, rotMatrix_z(1, self.phi)))
        self.rotMatrix_2 = np.dot(np.dot(X_90_2, rotMatrix_z(2, self.theta)), np.dot(X90_2, rotMatrix_z(2, self.phi)))
        self.rotMatrix_3 = np.dot(np.dot(X_90_3, rotMatrix_z(3, self.theta)), np.dot(X90_3, rotMatrix_z(3, self.phi)))

        #Rz = np.array([[np.cos(self.phi), np.sin(self.phi), 0], [-np.sin(self.phi), np.cos(self.phi), 0], [0, 0, 1]])
        #Rx = np.array([[1, 0, 0], [0, np.cos(self.theta), np.sin(self.theta)], [0, -np.sin(self.theta), np.cos(self.theta)]])
        #rotMatrix = np.dot(Rz, Rx)
        mat = np.array([self.x, self.y, self.z])
        #mat = np.dot(rotMatrix, mat)
        self.x = mat[0]
        self.y = mat[1]
        self.z = mat[2]
        # ShimTypes
        self.Y = np.sqrt(3/(4*np.pi)) * self.y
        self.X = np.sqrt(3/(4*np.pi)) * self.x
        self.Z = np.sqrt(3/(4*np.pi)) * self.z
        
        self.XY = 0.5*np.sqrt(15/np.pi) * self.x * self.y
        self.YZ = 0.5*np.sqrt(15/np.pi) * self.y * self.z
        self.Z2 = 0.25*np.sqrt(5/np.pi) * (2*self.z**2 - self.x**2 - self.y**2)
        self.XZ = 0.5*np.sqrt(15/np.pi) * self.x * self.z
        self.X2_Y2 = 0.25*np.sqrt(15/np.pi) * (self.x**2 - self.y**2)

        self.Y3 = 0.25*np.sqrt(35/(2*np.pi)) * (3*self.x**2 * self.y - self.y**3)
        self.XYZ = 0.5*np.sqrt(105/np.pi) * self.z * self.x * self.y
        self.YZ2 = 0.25*np.sqrt(21/(2*np.pi)) * self.y * (4*self.z**2 - self.x**2 - self.y**2)
        self.Z3 = 0.25*np.sqrt(7/np.pi)*self.z*(2*self.z**2 - 3*self.x**2 - 3*self.y**2)
        self.XZ2 = 0.25*np.sqrt(21/(2*np.pi)) * self.x * (4*self.z**2 - self.x**2 - self.y**2)
        self.ZX2_ZY2 = 0.25*np.sqrt(105/np.pi) * self.z * (self.x**2 - self.y**2)
        self.X3 = 0.25*np.sqrt(35/(2*np.pi)) * (self.x**3 - 3*self.y**2 * self.x) 

        self.Mfield = np.zeros(self.x.shape)
        self.freq = np.linspace(-self.sw/2.0, self.sw/2.0, self.npoints)
        self.lb = np.exp(-((10*np.pi*np.arange(self.npoints*self.OVERSAMPLE)/self.sw)**2) / (4.0 * np.log(2)))/self.N
        self.lb[self.npoints-(self.npoints+1)//2:] = 0
        self.lb[0] = self.lb[0]/2.0
        self.lbScale = np.sum(np.abs(self.lb)) * self.N / (self.npoints*self.OVERSAMPLE*100.0)

    def simulate(self, y1=0, z1=0, x1=0, xy=0, yz=0, z2=0, xz=0, x2_y2=0, y3=0, xyz=0, yz2=0, z3=0, xz2=0, zx2_zy2=0, x3=0, spinning=False):
        y1 += self.y1Game
        z1 += self.z1Game
        x1 += self.x1Game

        xy += self.xyGame
        yz += self.yzGame
        z2 += self.z2Game
        xz += self.xzGame
        x2_y2 += self.x2_y2Game

        y3 += self.y3Game
        xyz += self.xyzGame
        yz2 += self.yz2Game
        z3 += self.z3Game
        xz2 += self.xz2Game
        zx2_zy2 += self.zx2_zy2Game
        x3 += self.x3Game

        [y1, z1, x1] = np.dot(self.rotMatrix_1, [y1, z1, x1])
        [xy, yz, z2, xz, x2_y2] = np.dot(self.rotMatrix_2, [xy, yz, z2, xz, x2_y2])
        [y3, xyz, yz2, z3, xz2, zx2_zy2, x3] = np.dot(self.rotMatrix_3, [y3, xyz, yz2, z3, xz2, zx2_zy2, x3])

        self.Mfield = z1*self.Z + z2*self.Z2 + z3*self.Z3
        if not spinning:
            self.Mfield += x1*self.X + y1*self.Y
            self.Mfield += xz*self.XZ + xy*self.XY + yz*self.YZ + x2_y2*self.X2_Y2
            self.Mfield += xz2*self.XZ2 + yz2*self.YZ2 + zx2_zy2*self.ZX2_ZY2 + xyz*self.XYZ + x3*self.X3 + y3*self.Y3
        
        self.spectrum, tmp = np.histogram(self.Mfield, self.npoints*self.OVERSAMPLE, (-self.sw/2.0, self.sw/2.0))
        self.spectrum = np.fft.ifft(self.spectrum)*self.lb
        self.fidSurface = (np.sum(np.abs(self.spectrum)))/ (self.lbScale)        
        self.spectrum = np.fft.fft(self.spectrum)[::self.OVERSAMPLE]
        self.spectrum = np.real(self.spectrum)
        r  = UnivariateSpline(self.freq, self.spectrum-np.max(self.spectrum)/2.0, s=0).roots()
        self.fwhm = np.abs(max(r)-min(r))

    def resetGame(self):
        self.y1Game = 0.0
        self.z1Game = 0.0
        self.x1Game = 0.0
        
        self.xyGame = 0.0
        self.yzGame = 0.0
        self.z2Game = 0.0
        self.xzGame = 0.0
        self.x2_y2Game = 0.0

        self.y3Game = 0.0
        self.xyzGame = 0.0
        self.yz2Game = 0.0
        self.z3Game = 0.0
        self.xz2Game = 0.0
        self.zx2_zy2Game = 0.0
        self.x3Game = 0.0

    def startGame(self, order=4, zonly=False):
        self.resetGame()
        if zonly:
            self.z1Game = np.random.randint(-NSTEPS, NSTEPS+1)*Z1LIM/NSTEPS
            self.z2Game = np.random.randint(-NSTEPS, NSTEPS+1)*Z2LIM/NSTEPS
            self.z3Game = np.random.randint(-NSTEPS, NSTEPS+1)*Z3LIM/NSTEPS
            return
        self.z1Game = np.random.randint(-NSTEPS, NSTEPS+1)*Z1LIM/NSTEPS
        self.x1Game = np.random.randint(-NSTEPS, NSTEPS+1)*X1LIM/NSTEPS
        self.y1Game = np.random.randint(-NSTEPS, NSTEPS+1)*Y1LIM/NSTEPS
        if order > 1:
            self.z2Game = np.random.randint(-NSTEPS, NSTEPS+1)*Z2LIM/NSTEPS
            self.xzGame = np.random.randint(-NSTEPS, NSTEPS+1)*XZLIM/NSTEPS
            self.yzGame = np.random.randint(-NSTEPS, NSTEPS+1)*YZLIM/NSTEPS
            self.x2_y2Game = np.random.randint(-NSTEPS, NSTEPS+1)*X2_Y2LIM/NSTEPS
            self.xyGame = np.random.randint(-NSTEPS, NSTEPS+1)*XYLIM/NSTEPS
        if order > 2:
            self.z3Game = np.random.randint(-NSTEPS, NSTEPS+1)*Z3LIM/NSTEPS
            self.xz2Game = np.random.randint(-NSTEPS, NSTEPS+1)*XZ2LIM/NSTEPS
            self.yz2Game = np.random.randint(-NSTEPS, NSTEPS+1)*YZ2LIM/NSTEPS
            self.zx2_zy2Game = np.random.randint(-NSTEPS, NSTEPS+1)*ZX2_ZY2LIM/NSTEPS
            self.xyzGame = np.random.randint(-NSTEPS, NSTEPS+1)*XYZLIM/NSTEPS
            self.x3Game = np.random.randint(-NSTEPS, NSTEPS+1)*X3LIM/NSTEPS
            self.y3Game = np.random.randint(-NSTEPS, NSTEPS+1)*Y3LIM/NSTEPS
        

class ShimFrame(QtWidgets.QScrollArea):

    def __init__(self, parent):
        super(ShimFrame, self).__init__(parent)
        self.father = parent
        content = QtWidgets.QWidget(self)
        self.grid = QtWidgets.QGridLayout(content)

        self.grid.addWidget(QtWidgets.QLabel("FID Surface [%]:"), 0, 0)
        self.fidSurfaceLine = QtWidgets.QLineEdit(self)
        self.fidSurfaceLine.setAlignment(QtCore.Qt.AlignHCenter)
        self.fidSurfaceLine.setText("0.0")
        self.fidSurfaceLine.setReadOnly(True)
        self.grid.addWidget(self.fidSurfaceLine, 1, 0)
        
        self.grid.addWidget(QtWidgets.QLabel("FWHM:"), 0, 1)
        self.fwhmLine = QtWidgets.QLineEdit(self)
        self.fwhmLine.setAlignment(QtCore.Qt.AlignHCenter)
        self.fwhmLine.setText("0.0")
        self.fwhmLine.setReadOnly(True)
        self.grid.addWidget(self.fwhmLine, 1, 1)
        
        self.z1Label = QtWidgets.QLabel()
        self.grid.addWidget(self.z1Label, 2, 0)
        self.z1 = QtWidgets.QLineEdit(self)
        self.z1.setAlignment(QtCore.Qt.AlignHCenter)
        self.z1.setText("0.0")
        self.z1.returnPressed.connect(lambda: self.insertVal(self.z1, self.z1Scale, NSTEPS/Z1LIM))
        self.grid.addWidget(self.z1, 3, 0)
        self.z1Scale = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.z1Scale.setRange(-NSTEPS, NSTEPS)
        self.z1Scale.valueChanged.connect(lambda inp: self.setVal(inp, self.z1, NSTEPS/Z1LIM))
        self.grid.addWidget(self.z1Scale, 4, 0)

        self.x1Label = QtWidgets.QLabel()
        self.grid.addWidget(self.x1Label, 5, 0)
        self.x1 = QtWidgets.QLineEdit(self)
        self.x1.setAlignment(QtCore.Qt.AlignHCenter)
        self.x1.setText("0.0")
        self.x1.returnPressed.connect(lambda: self.insertVal(self.x1, self.x1Scale, NSTEPS/X1LIM))
        self.grid.addWidget(self.x1, 6, 0)
        self.x1Scale = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.x1Scale.setRange(-NSTEPS, NSTEPS)
        self.x1Scale.valueChanged.connect(lambda inp: self.setVal(inp, self.x1, NSTEPS/X1LIM))
        self.grid.addWidget(self.x1Scale, 7, 0)

        self.y1Label = QtWidgets.QLabel()
        self.grid.addWidget(self.y1Label, 8, 0)
        self.y1 = QtWidgets.QLineEdit(self)
        self.y1.setAlignment(QtCore.Qt.AlignHCenter)
        self.y1.setText("0.0")
        self.y1.returnPressed.connect(lambda: self.insertVal(self.y1, self.y1Scale, NSTEPS/Y1LIM))
        self.grid.addWidget(self.y1, 9, 0)
        self.y1Scale = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.y1Scale.setRange(-NSTEPS, NSTEPS)
        self.y1Scale.valueChanged.connect(lambda inp: self.setVal(inp, self.y1, NSTEPS/Y1LIM))
        self.grid.addWidget(self.y1Scale, 10, 0)

        self.z2Label = QtWidgets.QLabel()
        self.grid.addWidget(self.z2Label, 2, 1)
        self.z2 = QtWidgets.QLineEdit(self)
        self.z2.setAlignment(QtCore.Qt.AlignHCenter)
        self.z2.setText("0.0")
        self.z2.returnPressed.connect(lambda: self.insertVal(self.z2, self.z2Scale, NSTEPS/Z2LIM))
        self.grid.addWidget(self.z2, 3, 1)
        self.z2Scale = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.z2Scale.setRange(-NSTEPS, NSTEPS)
        self.z2Scale.valueChanged.connect(lambda inp: self.setVal(inp, self.z2, NSTEPS/Z2LIM))
        self.grid.addWidget(self.z2Scale, 4, 1)
        
        self.xzLabel = QtWidgets.QLabel()
        self.grid.addWidget(self.xzLabel, 5, 1)
        self.xz = QtWidgets.QLineEdit(self)
        self.xz.setAlignment(QtCore.Qt.AlignHCenter)
        self.xz.setText("0.0")
        self.xz.returnPressed.connect(lambda: self.insertVal(self.xz, self.xzScale, NSTEPS/XZLIM))
        self.grid.addWidget(self.xz, 6, 1)
        self.xzScale = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.xzScale.setRange(-NSTEPS, NSTEPS)
        self.xzScale.valueChanged.connect(lambda inp: self.setVal(inp, self.xz, NSTEPS/XZLIM))
        self.grid.addWidget(self.xzScale, 7, 1)

        self.yzLabel = QtWidgets.QLabel()
        self.grid.addWidget(self.yzLabel, 8, 1)
        self.yz = QtWidgets.QLineEdit(self)
        self.yz.setAlignment(QtCore.Qt.AlignHCenter)
        self.yz.setText("0.0")
        self.yz.returnPressed.connect(lambda: self.insertVal(self.yz, self.yzScale, NSTEPS/YZLIM))
        self.grid.addWidget(self.yz, 9, 1)
        self.yzScale = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.yzScale.setRange(-NSTEPS, NSTEPS)
        self.yzScale.valueChanged.connect(lambda inp: self.setVal(inp, self.yz, NSTEPS/YZLIM))
        self.grid.addWidget(self.yzScale, 10, 1)

        self.x2_y2Label = QtWidgets.QLabel()
        self.grid.addWidget(self.x2_y2Label, 11, 1)
        self.x2_y2 = QtWidgets.QLineEdit(self)
        self.x2_y2.setAlignment(QtCore.Qt.AlignHCenter)
        self.x2_y2.setText("0.0")
        self.x2_y2.returnPressed.connect(lambda: self.insertVal(self.x2_y2, self.x2_y2Scale, NSTEPS/X2_Y2LIM))
        self.grid.addWidget(self.x2_y2, 12, 1)
        self.x2_y2Scale = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.x2_y2Scale.setRange(-NSTEPS, NSTEPS)
        self.x2_y2Scale.valueChanged.connect(lambda inp: self.setVal(inp, self.x2_y2, NSTEPS/X2_Y2LIM))
        self.grid.addWidget(self.x2_y2Scale, 13, 1)

        self.xyLabel = QtWidgets.QLabel()
        self.grid.addWidget(self.xyLabel, 14, 1)
        self.xy = QtWidgets.QLineEdit(self)
        self.xy.setAlignment(QtCore.Qt.AlignHCenter)
        self.xy.setText("0.0")
        self.xy.returnPressed.connect(lambda: self.insertVal(self.xy, self.xyScale, NSTEPS/XYLIM))
        self.grid.addWidget(self.xy, 15, 1)
        self.xyScale = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.xyScale.setRange(-NSTEPS, NSTEPS)
        self.xyScale.valueChanged.connect(lambda inp: self.setVal(inp, self.xy, NSTEPS/XYLIM))
        self.grid.addWidget(self.xyScale, 16, 1)
        
        self.z3Label = QtWidgets.QLabel()
        self.grid.addWidget(self.z3Label, 2, 2)
        self.z3 = QtWidgets.QLineEdit(self)
        self.z3.setAlignment(QtCore.Qt.AlignHCenter)
        self.z3.setText("0.0")
        self.z3.returnPressed.connect(lambda: self.insertVal(self.z3, self.z3Scale, NSTEPS/Z3LIM))
        self.grid.addWidget(self.z3, 3, 2)
        self.z3Scale = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.z3Scale.setRange(-NSTEPS, NSTEPS)
        self.z3Scale.valueChanged.connect(lambda inp: self.setVal(inp, self.z3, NSTEPS/Z3LIM))
        self.grid.addWidget(self.z3Scale, 4, 2)

        self.xz2Label = QtWidgets.QLabel()
        self.grid.addWidget(self.xz2Label, 5, 2)
        self.xz2 = QtWidgets.QLineEdit(self)
        self.xz2.setAlignment(QtCore.Qt.AlignHCenter)
        self.xz2.setText("0.0")
        self.xz2.returnPressed.connect(lambda: self.insertVal(self.xz2, self.xz2Scale, NSTEPS/XZ2LIM))
        self.grid.addWidget(self.xz2, 6, 2)
        self.xz2Scale = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.xz2Scale.setRange(-NSTEPS, NSTEPS)
        self.xz2Scale.valueChanged.connect(lambda inp: self.setVal(inp, self.xz2, NSTEPS/XZ2LIM))
        self.grid.addWidget(self.xz2Scale, 7, 2)

        self.yz2Label = QtWidgets.QLabel()
        self.grid.addWidget(self.yz2Label, 8, 2)
        self.yz2 = QtWidgets.QLineEdit(self)
        self.yz2.setAlignment(QtCore.Qt.AlignHCenter)
        self.yz2.setText("0.0")
        self.yz2.returnPressed.connect(lambda: self.insertVal(self.yz2, self.yz2Scale, NSTEPS/YZ2LIM))
        self.grid.addWidget(self.yz2, 9, 2)
        self.yz2Scale = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.yz2Scale.setRange(-NSTEPS, NSTEPS)
        self.yz2Scale.valueChanged.connect(lambda inp: self.setVal(inp, self.yz2, NSTEPS/YZ2LIM))
        self.grid.addWidget(self.yz2Scale, 10, 2)
        
        self.zx2_zy2Label = QtWidgets.QLabel()
        self.grid.addWidget(self.zx2_zy2Label, 11, 2)
        self.zx2_zy2 = QtWidgets.QLineEdit(self)
        self.zx2_zy2.setAlignment(QtCore.Qt.AlignHCenter)
        self.zx2_zy2.setText("0.0")
        self.zx2_zy2.returnPressed.connect(lambda: self.insertVal(self.zx2_zy2, self.zx2_zy2Scale, NSTEPS/ZX2_ZY2LIM))
        self.grid.addWidget(self.zx2_zy2, 12, 2)
        self.zx2_zy2Scale = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.zx2_zy2Scale.setRange(-NSTEPS, NSTEPS)
        self.zx2_zy2Scale.valueChanged.connect(lambda inp: self.setVal(inp, self.zx2_zy2, NSTEPS/ZX2_ZY2LIM))
        self.grid.addWidget(self.zx2_zy2Scale, 13, 2)
         
        self.xyzLabel = QtWidgets.QLabel()
        self.grid.addWidget(self.xyzLabel, 14, 2)
        self.xyz = QtWidgets.QLineEdit(self)
        self.xyz.setAlignment(QtCore.Qt.AlignHCenter)
        self.xyz.setText("0.0")
        self.xyz.returnPressed.connect(lambda: self.insertVal(self.xyz, self.xyzScale, NSTEPS/XYZLIM))
        self.grid.addWidget(self.xyz, 15, 2)
        self.xyzScale = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.xyzScale.setRange(-NSTEPS, NSTEPS)
        self.xyzScale.valueChanged.connect(lambda inp: self.setVal(inp, self.xyz, NSTEPS/XYZLIM))
        self.grid.addWidget(self.xyzScale, 16, 2)

        self.x3Label = QtWidgets.QLabel()
        self.grid.addWidget(self.x3Label, 17, 2)
        self.x3 = QtWidgets.QLineEdit(self)
        self.x3.setAlignment(QtCore.Qt.AlignHCenter)
        self.x3.setText("0.0")
        self.x3.returnPressed.connect(lambda: self.insertVal(self.x3, self.x3Scale, NSTEPS/X3LIM))
        self.grid.addWidget(self.x3, 18, 2)
        self.x3Scale = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.x3Scale.setRange(-NSTEPS, NSTEPS)
        self.x3Scale.valueChanged.connect(lambda inp: self.setVal(inp, self.x3, NSTEPS/X3LIM))
        self.grid.addWidget(self.x3Scale, 19, 2)

        self.y3Label = QtWidgets.QLabel()
        self.grid.addWidget(self.y3Label, 20, 2)
        self.y3 = QtWidgets.QLineEdit(self)
        self.y3.setAlignment(QtCore.Qt.AlignHCenter)
        self.y3.setText("0.0")
        self.y3.returnPressed.connect(lambda: self.insertVal(self.y3, self.y3Scale, NSTEPS/Y3LIM))
        self.grid.addWidget(self.y3, 21, 2)
        self.y3Scale = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.y3Scale.setRange(-NSTEPS, NSTEPS)
        self.y3Scale.valueChanged.connect(lambda inp: self.setVal(inp, self.y3, NSTEPS/Y3LIM))
        self.grid.addWidget(self.y3Scale, 22, 2)
        
        #Make widget lists
         
        self.lineEditsFirst = [self.y1 , self.x1 , self.z1]
        self.lineEditsSecond = [self.z2 , self.xz, self.yz , self.xy, self.x2_y2]
        self.lineEditsThird = [self.y3 , self.xyz , self.yz2 , self.z3 ,
                          self.xz2 , self.xz2 , self.zx2_zy2 , self.x3]
        self.lineEditsZ = [self.z1 , self.z2, self.z3]
        self.lineEditsNonZ = [self.y1 , self.x1 , self.xz, 
                         self.yz , self.xy, self.x2_y2,
                         self.y3 , self.xyz , self.yz2 ,
                          self.xz2 , self.xz2 , self.zx2_zy2 , self.x3]
                          
        self.sliderFirst = [self.y1Scale , self.x1Scale , self.z1Scale]
        self.sliderSecond = [self.z2Scale , self.xzScale, self.yzScale , self.xyScale, self.x2_y2Scale]
        self.sliderThird = [self.y3Scale , self.xyzScale , self.yz2Scale , self.z3Scale ,
                          self.xz2Scale , self.xz2Scale , self.zx2_zy2Scale , self.x3Scale]

        self.sliderZ = [self.z1Scale , self.z2Scale, self.z3Scale]
        self.sliderNonZ = [self.y1Scale , self.x1Scale , self.xzScale, 
                         self.yzScale , self.xyScale, self.x2_y2Scale,
                         self.y3Scale , self.xyzScale , self.yz2Scale ,
                          self.xz2Scale , self.xz2Scale , self.zx2_zy2Scale , self.x3Scale]
        
        
        
        QtCore.QTimer.singleShot(100, self.resizeAll)
        self.setTitleText()

    def setTitleText(self, y1='', z1='', x1='', xy='', yz='', z2='', xz= '', x2_y2='', y3='', xyz='', yz2='', z3='', xz2='', zx2_zy2='', x3=''):
        self.y1Label.setText('Y: ' + str(y1))
        self.z1Label.setText('Z: ' + str(z1))
        self.x1Label.setText('X: ' + str(x1))
        self.xyLabel.setText('XY: ' + str(xy))
        self.yzLabel.setText('YZ: ' + str(yz))
        self.z2Label.setText('Z2: ' + str(z2))
        self.xzLabel.setText('XZ: ' + str(xz))
        self.x2_y2Label.setText('X2-Y2: ' + str(x2_y2))
        self.y3Label.setText('Y3: ' + str(y3))
        self.xyzLabel.setText('XYZ: ' + str(xyz))
        self.yz2Label.setText('YZ2: ' + str(yz2))
        self.z3Label.setText('Z3: ' + str(z3))
        self.xz2Label.setText('XZ2: ' + str(xz2))
        self.zx2_zy2Label.setText('(X2-Y2)Z: ' + str(zx2_zy2))
        self.x3Label.setText('X3: ' + str(x3))
        
    def resizeAll(self):
        self.setMinimumWidth(self.grid.sizeHint().width() + self.verticalScrollBar().sizeHint().width())
        
    def insertVal(self, func, scaleFunc, scale):
        inp = safeEval(func.text())
        scaleFunc.blockSignals(True)
        if inp is None:
            scaleFunc.setSliderPosition(0)
        else:
            scaleFunc.setSliderPosition(int(round(inp*scale)))
        scaleFunc.blockSignals(False)
        self.sim()
            
    def setVal(self, val, func, scale):
        func.setText(str(val/scale))
        self.sim()
        
    def sim(self, *args):
        y1 = safeEval(self.y1.text())
        z1 = safeEval(self.z1.text())
        x1 = safeEval(self.x1.text())
        xy = safeEval(self.xy.text())
        yz = safeEval(self.yz.text())
        z2 = safeEval(self.z2.text())
        xz = safeEval(self.xz.text())
        x2_y2 = safeEval(self.x2_y2.text())
        y3 = safeEval(self.y3.text())
        xyz = safeEval(self.xyz.text())
        yz2 = safeEval(self.yz2.text())
        z3 = safeEval(self.z3.text())
        xz2 = safeEval(self.xz2.text())
        zx2_zy2 = safeEval(self.zx2_zy2.text())
        x3 = safeEval(self.x3.text())
        self.father.sim(z1=z1, x1=x1, y1=y1, xy=xy, yz=yz, z2=z2, xz=xz, x2_y2=x2_y2, y3=y3, xyz=xyz, yz2=yz2, z3=z3, xz2=xz2, zx2_zy2=zx2_zy2, x3=x3)

    def resetShims(self):
        lineedits = [self.y1, self.z1, self.x1, self.xz, self.yz, self.z2, self.xz, self.x2_y2, self.y3, self.xyz, self.yz2, self.z3, self.xz2, self.zx2_zy2, self.x3]
        scale = [self.y1Scale, self.z1Scale, self.x1Scale, self.xyScale, self.yzScale, self.z2Scale, self.xzScale, self.x2_y2Scale, self.y3Scale, self.xyzScale, self.yz2Scale, self.z3Scale, self.xz2Scale, self.zx2_zy2Scale, self.x3Scale]
        for i in range(len(lineedits)):
            lineedits[i].setText("0.0")
            scale[i].blockSignals(True)
            scale[i].setValue(0)
            scale[i].blockSignals(False)
        self.sim()

    def setValues(self, fidSurface, fwhm):
        if fidSurface > 95:
            self.fidSurfaceLine.setStyleSheet("color: rgb(0, 137, 0);")
        else:
            self.fidSurfaceLine.setStyleSheet("color: rgb(0, 0, 0);")
        self.fidSurfaceLine.setText("%.2f" % fidSurface)
        self.fwhmLine.setText("%.2f" % fwhm)


class SettingsFrame(QtWidgets.QWidget):

    def __init__(self, parent, shimSim):
        super(SettingsFrame, self).__init__(parent)
        self.father = parent
        self.shimSim = shimSim
        grid = QtWidgets.QGridLayout(self)
        grid.addWidget(QtWidgets.QLabel("Settings:"), 0, 0)
        grid.addWidget(QtWidgets.QLabel("Diameter [mm]:"), 1, 0)
        self.diameter = QtWidgets.QLineEdit(self)
        self.diameter.setAlignment(QtCore.Qt.AlignHCenter)
        self.diameter.setText(str(self.shimSim.r*2))
        self.diameter.returnPressed.connect(self.setDimensions)
        grid.addWidget(self.diameter, 1, 1)
        grid.addWidget(QtWidgets.QLabel("Length [mm]"), 2, 0)
        self.length = QtWidgets.QLineEdit(self)
        self.length.setAlignment(QtCore.Qt.AlignHCenter)
        self.length.setText(str(self.shimSim.l))
        self.length.returnPressed.connect(self.setDimensions)
        grid.addWidget(self.length, 2, 1)
        self.spinningButton = QtWidgets.QCheckBox('Spinning', self)
        self.spinningButton.stateChanged.connect(self.father.shimFrame.sim)
        grid.addWidget(self.spinningButton, 3, 0)

        
        grid.addWidget(QtWidgets.QLabel("Rotation:"), 0, 2)
        grid.addWidget(QtWidgets.QLabel("Theta [degrees]:"), 1, 2)
        self.theta = QtWidgets.QLineEdit(self)
        self.theta.setAlignment(QtCore.Qt.AlignHCenter)
        self.theta.setText(str(self.shimSim.theta/np.pi*180.0))
        self.theta.returnPressed.connect(self.setAngles)
        grid.addWidget(self.theta, 1, 3)
        grid.addWidget(QtWidgets.QLabel("Phi [degrees]:"), 2, 2)
        self.phi = QtWidgets.QLineEdit(self)
        self.phi.setAlignment(QtCore.Qt.AlignHCenter)
        self.phi.setText(str(self.shimSim.phi/np.pi*180.0))
        self.phi.returnPressed.connect(self.setAngles)
        grid.addWidget(self.phi, 2, 3)
        
        grid.addWidget(QtWidgets.QLabel("# points:"), 1, 4)
        self.nsamples = QtWidgets.QLineEdit(self)
        self.nsamples.setAlignment(QtCore.Qt.AlignHCenter)
        self.nsamples.setText(str(self.shimSim.N))
        self.nsamples.returnPressed.connect(self.setPoints)
        grid.addWidget(self.nsamples, 1, 5)

        grid.addWidget(QtWidgets.QLabel("Sweepwidth:"), 1, 6)
        self.sw = QtWidgets.QLineEdit(self)
        self.sw.setAlignment(QtCore.Qt.AlignHCenter)
        self.sw.setText(str(self.shimSim.sw))
        self.sw.returnPressed.connect(self.setSw)
        grid.addWidget(self.sw, 1, 7)
        
        grid.addWidget(QtWidgets.QLabel("Game:"), 0, 8)
        self.game1Button = QtWidgets.QPushButton('1st order', self)
        self.game1Button.clicked.connect(lambda: self.father.startGame(1))
        grid.addWidget(self.game1Button, 1, 8)
        self.game2Button = QtWidgets.QPushButton('2nd order', self)
        self.game2Button.clicked.connect(lambda: self.father.startGame(2))
        grid.addWidget(self.game2Button, 2, 8)
        self.game3Button = QtWidgets.QPushButton('3rd order', self)
        self.game3Button.clicked.connect(lambda: self.father.startGame(3))
        grid.addWidget(self.game3Button, 3, 8)
        self.gameZButton = QtWidgets.QPushButton('Z* only', self)
        self.gameZButton.clicked.connect(lambda: self.father.startGame(zonly=True))
        grid.addWidget(self.gameZButton, 4, 8)
        self.gameStopButton = QtWidgets.QPushButton('Stop', self)
        self.gameStopButton.clicked.connect(self.father.resetGame)
        grid.addWidget(self.gameStopButton, 5, 8)
        
        self.resultsButton = QtWidgets.QCheckBox('Results', self)
        self.resultsButton.stateChanged.connect(self.results)
        grid.addWidget(self.resultsButton, 0, 9)
        
        self.resetButton = QtWidgets.QPushButton('Reset Shims', self)
        self.resetButton.clicked.connect(self.father.resetShims)
        grid.addWidget(self.resetButton, 1, 9)
        
        grid.setColumnStretch(10, 1)
        grid.setRowStretch(10, 1)

    def setDimensions(self):
        diameter = safeEval(self.diameter.text())
        length = safeEval(self.length.text())
        if diameter is not None:
            self.shimSim.r = diameter/2.0
        if length is not None:
            self.shimSim.l = length
        self.shimSim.setupGrid()
        self.father.shimFrame.sim()

    def setAngles(self):
        theta = safeEval(self.theta.text())
        phi = safeEval(self.phi.text())
        if theta is not None:
            self.shimSim.theta = theta/180.0*np.pi
        if phi is not None:
            self.shimSim.phi = phi/180.0*np.pi
        self.shimSim.setupGrid()
        self.father.shimFrame.sim()

    def setPoints(self):
        nsamples = safeEval(self.nsamples.text())
        if nsamples is not None:
            self.shimSim.N = int(round(nsamples))
        self.shimSim.setupGrid()
        self.father.shimFrame.sim()

    def setSw(self):
        sw = safeEval(self.sw.text())
        if sw is not None:
            self.shimSim.sw = int(round(sw))
        self.shimSim.setupGrid()
        self.father.shimFrame.sim()
        
    def results(self, val):
        if val:
            self.father.shimFrame.setTitleText(-self.shimSim.y1Game, -self.shimSim.z1Game, -self.shimSim.x1Game, -self.shimSim.xyGame, -self.shimSim.yzGame, -self.shimSim.z2Game, -self.shimSim.xzGame, -self.shimSim.x2_y2Game, -self.shimSim.y3Game, -self.shimSim.xyzGame, -self.shimSim.yz2Game, -self.shimSim.z3Game, -self.shimSim.xz2Game, -self.shimSim.zx2_zy2Game, -self.shimSim.x3Game)
        else:
            self.father.shimFrame.setTitleText()


class ShimPlotFrame(Plot1DFrame):

    def __init__(self, root, fig, canvas):
        super(ShimPlotFrame, self).__init__(root, fig, canvas)
        self.canvas.mpl_connect('button_press_event', self.buttonPress)
        self.canvas.mpl_connect('button_release_event', self.buttonRelease)
        self.canvas.mpl_connect('motion_notify_event', self.pan)
        self.canvas.mpl_connect('scroll_event', self.scroll)
        self.canvas.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.canvas.setFocus()
        self.xmaxlim= None
        self.xminlim= None
        self.ymaxlim= None
        self.yminlim= None

    def setData(self, xdata, ydata):
        self.xdata = xdata
        self.ydata = ydata

    def plotReset(self, xReset=True, yReset=True):  # set the plot limits to min and max values
        miny = min(self.ydata)
        maxy = max(self.ydata)
        differ = 0.05 * (maxy - miny)  # amount to add to show all datapoints (10%)
        if yReset:
            self.yminlim = miny - differ
            self.ymaxlim = maxy + differ
        axMult = 1.0 
        if xReset:
            self.xminlim = min(self.xdata * axMult)
            self.xmaxlim = max(self.xdata * axMult)
        self.ax.set_xlim(self.xmaxlim, self.xminlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        
    def showFid(self):
        self.ax.cla()
        self.ax.plot(self.xdata, self.ydata)
        if self.xmaxlim is None:
            self.plotReset()
        self.ax.set_xlim(self.xmaxlim, self.xminlim)
        self.ax.set_ylim(self.yminlim, self.ymaxlim)
        self.canvas.draw()        

class MainProgram(QtWidgets.QMainWindow):

    def __init__(self, root):
        super(MainProgram, self).__init__()
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.main_widget = QtWidgets.QWidget(self)
        self.mainFrame = QtWidgets.QGridLayout(self.main_widget)
        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)
        self.menubar = self.menuBar()
        self.filemenu = QtWidgets.QMenu('&File', self)
        self.menubar.addMenu(self.filemenu)
        self.savefigAct = self.filemenu.addAction('Export Figure', self.saveFigure, QtGui.QKeySequence.Print)
        self.savefigAct.setToolTip('Export as Figure')
        self.savedatAct = self.filemenu.addAction('Export Data', self.saveData)
        self.savedatAct.setToolTip('Export as text')
        self.quitAct = self.filemenu.addAction('&Quit', self.fileQuit, QtGui.QKeySequence.Quit)
        self.quitAct.setToolTip('Close Shimpanzee')

        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        self.ax = self.fig.gca()
        self.mainFrame.addWidget(self.canvas, 0, 0)
        self.mainFrame.setColumnStretch(0, 1)
        self.mainFrame.setRowStretch(0, 1)
        self.shimFrame = ShimFrame(self)
        self.mainFrame.addWidget(self.shimFrame, 0, 1, 2, 1)
        self.shimSim = ShimSim()
        self.settingsFrame = SettingsFrame(self, self.shimSim)
        self.mainFrame.addWidget(self.settingsFrame, 1, 0)
        self.shimPlotFrame = ShimPlotFrame(self, self.fig, self.canvas)
        self.sim()
        self.shimPlotFrame.plotReset()

    def sim(self, *args, **kwargs):
        self.shimSim.simulate(*args, spinning=self.settingsFrame.spinningButton.isChecked(), **kwargs)
        self.shimPlotFrame.setData(self.shimSim.freq, self.shimSim.spectrum)
        self.shimPlotFrame.showFid()
        self.shimFrame.setValues(self.shimSim.fidSurface, self.shimSim.fwhm)

    def startGame(self, *args, **kwargs):
        self.shimSim.startGame(*args, **kwargs)
       
        
        
        if len(args) > 0:   
            if args[0] == 1:
                for widget in self.shimFrame.lineEditsFirst + self.shimFrame.sliderFirst:
                    widget.setDisabled(False)
                for widget in self.shimFrame.lineEditsSecond + self.shimFrame.lineEditsThird + self.shimFrame.sliderSecond + self.shimFrame.sliderThird:
                    widget.setDisabled(True)
            elif args[0] == 2:
                for widget in self.shimFrame.lineEditsFirst + self.shimFrame.lineEditsSecond + self.shimFrame.sliderFirst + self.shimFrame.sliderSecond:
                    widget.setDisabled(False)
                for widget in self.shimFrame.lineEditsThird + self.shimFrame.sliderThird:
                    widget.setDisabled(True)
            elif args[0] == 3:
                for widget in self.shimFrame.lineEditsFirst + self.shimFrame.lineEditsSecond + self.shimFrame.lineEditsThird + self.shimFrame.sliderFirst + self.shimFrame.sliderSecond + self.shimFrame.sliderThird:
                    widget.setDisabled(False)
        elif kwargs['zonly']:
            for widget in self.shimFrame.lineEditsZ + self.shimFrame.sliderZ:
                widget.setDisabled(False)
            for widget in self.shimFrame.lineEditsNonZ + self.shimFrame.sliderNonZ:
                widget.setDisabled(True)

        
        
        
        
        self.resetShims()
        self.shimFrame.sim()

    def resetGame(self):
        self.shimSim.resetGame()
        for widget in self.shimFrame.lineEditsFirst + self.shimFrame.lineEditsSecond + self.shimFrame.lineEditsThird + self.shimFrame.sliderFirst + self.shimFrame.sliderSecond + self.shimFrame.sliderThird:
            widget.setDisabled(False)
        self.shimFrame.sim()

    def resetShims(self):
        self.shimFrame.resetShims()

    def saveFigure(self):
        f = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', 'Spectrum.png' ,filter = '(*.png)')
        if type(f) is tuple:
            f = f[0]        
        if f:
            dpi = 150
            self.fig.savefig(f, format='png', dpi=dpi)

    def saveData(self):
        f = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', 'Spectrum.txt' ,filter = '(*.txt)')
        if type(f) is tuple:
            f = f[0]        
        if f:
            data = np.zeros((len(self.shimSim.spectrum),2))
            data[:,0] = self.shimPlotFrame.xdata
            data[:,1] =  self.shimSim.spectrum
            np.savetxt(f,data)

    def fileQuit(self):
        self.close()
 
if __name__ == '__main__':
    mainProgram = MainProgram(root)
    mainProgram.setWindowTitle("Shimpanzee")
    mainProgram.show()
    splash.finish(mainProgram)
    sys.exit(root.exec_())
    

