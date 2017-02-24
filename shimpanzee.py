import sip
import sys
sip.setapi('QString', 2)
try:
    from PyQt4 import QtGui, QtCore
    from PyQt4 import QtGui as QtWidgets
    QT = 4
except ImportError:
    from PyQt5 import QtGui, QtCore, QtWidgets
    QT = 5
import matplotlib
if QT ==4:
    matplotlib.use('Qt4Agg')
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
else:
    matplotlib.use('Qt5Agg')
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from safeEval import safeEval
from spectrumFrame_Shim import Plot1DFrame
from scipy.interpolate import UnivariateSpline


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

Z4LIM = 0.25
Z5LIM = 0.05

NSTEPS = 1000


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
        Rz = np.array([[np.cos(self.phi), np.sin(self.phi), 0], [-np.sin(self.phi), np.cos(self.phi), 0], [0, 0, 1]])
        Rx = np.array([[1, 0, 0], [0, np.cos(self.theta), np.sin(self.theta)], [0, -np.sin(self.theta), np.cos(self.theta)]])
        rotMatrix = np.dot(Rz, Rx)
        mat = np.array([self.x, self.y, self.z])
        mat = np.dot(rotMatrix, mat)
        self.x = mat[0]
        self.y = mat[1]
        self.z = mat[2]
        # ShimTypes
        self.Z = self.z
        self.X = self.x
        self.Y = self.y
        
        self.Z2 = self.z**2 - 0.5*(self.x**2 + self.y**2)
        self.XZ = 3 * self.x * self.z
        self.YZ = 3 * self.y * self.z
        self.X2_Y2 = 3 * (self.x**2 - self.y**2)
        self.XY = 6 * self.x * self.y
        
        self.Z3 = self.z**3 - 3.0/2.0*(self.x**2 + self.y**2)*self.z
        self.XZ2 = 6 * self.x * (self.z**2 - (self.x**2 + self.y**2)/4.0)
        self.YZ2 = 6 * self.y * (self.z**2 - (self.x**2 + self.y**2)/4.0)
        self.ZX2_ZY2 = 15 * self.z * (self.x**2 - self.y**2)
        self.XYZ = 30 * self.z * self.x * self.y
        self.X3 = 15 * self.x**3 - 45 * self.x*self.y**2
        self.Y3 = 45 * self.x**2 * self.y - 15*self.y**3

        self.Z4 = self.z**4 - 3*self.z**2*(self.x**2 + self.y**2) + 3.0/8.0*(self.x**2 + self.y**2)**2
        self.Z5 = self.z**5 - 5*self.z**3*(self.x**2 + self.y**2) + 15.0/8.0*self.z*(self.x**2 + self.y**2)**2

        self.Mfield = np.zeros(self.x.shape)
        self.freq = np.linspace(-self.sw/2.0, self.sw/2.0, self.npoints)
        self.lb = np.exp(-((10*np.pi*np.arange(self.npoints*self.OVERSAMPLE)/self.sw)**2) / (4.0 * np.log(2)))
        self.lb[self.npoints-(self.npoints+1)//2:] = 0
        self.lb[0] = self.lb[0]/2.0

    def simulate(self, z1=0, x1=0, y1=0, z2=0, xz=0, yz=0, x2_y2=0, xy=0, z3=0, xz2=0, yz2=0, zx2_zy2=0, xyz=0, x3=0, y3=0, z4=0, z5=0):
        z1 += self.z1Game
        x1 += self.x1Game
        y1 += self.y1Game

        z2 += self.z2Game
        xz += self.xzGame
        yz += self.yzGame
        x2_y2 += self.x2_y2Game
        xy += self.xyGame

        z3 += self.z3Game
        xz2 += self.xz2Game
        yz2 += self.yz2Game
        zx2_zy2 += self.zx2_zy2Game
        xyz += self.xyzGame
        x3 += self.x3Game
        y3 += self.y3Game

        z4 += self.z4Game
        z5 += self.z5Game

        self.Mfield = x1*self.X + y1*self.Y + z1*self.Z
        self.Mfield += z2*self.Z2 + xz*self.XZ + xy*self.XY + yz*self.YZ + x2_y2*self.X2_Y2
        self.Mfield += z3*self.Z3 + xz2*self.XZ2 + yz2*self.YZ2 + zx2_zy2*self.ZX2_ZY2 + xyz*self.XYZ + x3*self.X3 + y3*self.Y3
        self.Mfield += z4*self.Z4 + z5*self.Z5  
        
        self.spectrum, tmp = np.histogram(self.Mfield, self.npoints*self.OVERSAMPLE, (-self.sw/2.0, self.sw/2.0), density=True)
        self.spectrum = np.fft.fft(np.fft.ifft(self.spectrum)*self.lb)
        self.spectrum = self.spectrum[::self.OVERSAMPLE]
        self.fidSurface = (np.sum(np.abs(self.spectrum))+np.abs(self.spectrum[0])*0.5)*self.sw
        self.spectrum = np.real(self.spectrum)
        r  = UnivariateSpline(self.freq, self.spectrum-np.max(self.spectrum)/2.0, s=0).roots()
        self.fwhm = np.abs(max(r)-min(r))

    def resetGame(self):
        self.z1Game = 0.0
        self.x1Game = 0.0
        self.y1Game = 0.0
        
        self.z2Game = 0.0
        self.xzGame = 0.0
        self.yzGame = 0.0
        self.x2_y2Game = 0.0
        self.xyGame = 0.0

        self.z3Game = 0.0
        self.xz2Game = 0.0
        self.yz2Game = 0.0
        self.zx2_zy2Game = 0.0
        self.xyzGame = 0.0
        self.x3Game = 0.0
        self.y3Game = 0.0

        self.z4Game = 0.0
        self.z5Game = 0.0

    def startGame(self, order=4, zonly=False):
        self.resetGame()
        if zonly:
            self.z1Game = np.random.randint(-NSTEPS, NSTEPS+1)*Z1LIM/NSTEPS
            self.z2Game = np.random.randint(-NSTEPS, NSTEPS+1)*Z2LIM/NSTEPS
            self.z3Game = np.random.randint(-NSTEPS, NSTEPS+1)*Z3LIM/NSTEPS
            self.z4Game = np.random.randint(-NSTEPS, NSTEPS+1)*Z4LIM/NSTEPS
            self.z5Game = np.random.randint(-NSTEPS, NSTEPS+1)*Z5LIM/NSTEPS
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
        if order > 3:
            self.z4Game = np.random.randint(-NSTEPS, NSTEPS+1)*Z4LIM/NSTEPS
            self.z5Game = np.random.randint(-NSTEPS, NSTEPS+1)*Z5LIM/NSTEPS
        

class ShimFrame(QtWidgets.QScrollArea):

    def __init__(self, parent):
        super(ShimFrame, self).__init__(parent)
        self.father = parent
        content = QtWidgets.QWidget(self)
        self.grid = QtWidgets.QGridLayout(content)

        self.grid.addWidget(QtWidgets.QLabel("FID Surface:"), 0, 0)
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
        
        self.z4Label = QtWidgets.QLabel()
        self.grid.addWidget(self.z4Label, 2, 3)
        self.z4 = QtWidgets.QLineEdit(self)
        self.z4.setAlignment(QtCore.Qt.AlignHCenter)
        self.z4.setText("0.0")
        self.z4.returnPressed.connect(lambda: self.insertVal(self.z4, self.z4Scale, NSTEPS/Z4LIM))
        self.grid.addWidget(self.z4, 3, 3)
        self.z4Scale = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.z4Scale.setRange(-NSTEPS, NSTEPS)
        self.z4Scale.valueChanged.connect(lambda inp: self.setVal(inp, self.z4, NSTEPS/Z4LIM))
        self.grid.addWidget(self.z4Scale, 4, 3)

        self.z5Label = QtWidgets.QLabel()
        self.grid.addWidget(self.z5Label, 5, 3)
        self.z5 = QtWidgets.QLineEdit(self)
        self.z5.setAlignment(QtCore.Qt.AlignHCenter)
        self.z5.setText("0.0")
        self.z5.returnPressed.connect(lambda: self.insertVal(self.z5, self.z5Scale, NSTEPS/Z5LIM))
        self.grid.addWidget(self.z5, 6, 3)
        self.z5Scale = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.z5Scale.setRange(-NSTEPS, NSTEPS)
        self.z5Scale.valueChanged.connect(lambda inp: self.setVal(inp, self.z5, NSTEPS/Z5LIM))
        self.grid.addWidget(self.z5Scale, 7, 3)
        QtCore.QTimer.singleShot(100, self.resizeAll)
        self.setTitleText()

    def setTitleText(self, z1='', x1='', y1='', z2='', xz='', yz='', x2_y2='', xy='', z3='', xz2='', yz2='', zx2_zy2='', xyz='', x3='', y3='', z4='', z5=''):
        self.z1Label.setText('Z: ' + str(z1))
        self.x1Label.setText('X: ' + str(x1))
        self.y1Label.setText('Y: ' + str(y1))
        self.z2Label.setText('Z2: ' + str(z2))
        self.xzLabel.setText('XZ: ' + str(xz))
        self.yzLabel.setText('YZ: ' + str(yz))
        self.x2_y2Label.setText('X2-Y2: ' + str(x2_y2))
        self.xyLabel.setText('XY: ' + str(xy))
        self.z3Label.setText('Z3: ' + str(z3))
        self.xz2Label.setText('XZ2: ' + str(xz2))
        self.yz2Label.setText('YZ2: ' + str(yz2))
        self.zx2_zy2Label.setText('(X2-Y2)Z: ' + str(zx2_zy2))
        self.xyzLabel.setText('XYZ: ' + str(xyz))
        self.x3Label.setText('X3: ' + str(x3))
        self.y3Label.setText('Y3: ' + str(y3))
        self.z4Label.setText('Z4: ' + str(z4))
        self.z5Label.setText('Z5: ' + str(z5))
        
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
        z1 = safeEval(self.z1.text())
        x1 = safeEval(self.x1.text())
        y1 = safeEval(self.y1.text())
        z2 = safeEval(self.z2.text())
        xz = safeEval(self.xz.text())
        yz = safeEval(self.yz.text())
        x2_y2 = safeEval(self.x2_y2.text())
        xy = safeEval(self.xy.text())
        z3 = safeEval(self.z3.text())
        xz2 = safeEval(self.xz2.text())
        yz2 = safeEval(self.yz2.text())
        zx2_zy2 = safeEval(self.zx2_zy2.text())
        xyz = safeEval(self.xyz.text())
        x3 = safeEval(self.x3.text())
        y3 = safeEval(self.y3.text())
        z4 = safeEval(self.z4.text())
        z5 = safeEval(self.z5.text())
        self.father.sim(z1=z1, x1=x1, y1=y1, z2=z2, xz=xz, yz=yz, x2_y2=x2_y2, xy=xy, z3=z3, xz2=xz2, yz2=yz2, zx2_zy2=zx2_zy2, xyz=xyz, x3=x3, y3=y3, z4=z4, z5=z5)

    def resetShims(self):
        lineedits = [self.z1, self.x1, self.y1, self.z2, self.xz, self.yz, self.x2_y2, self.xy, self.z3, self.xz2, self.yz2, self.zx2_zy2, self.xyz, self.x3, self.y3, self.z4, self.z5]
        scale = [self.z1Scale, self.x1Scale, self.y1Scale, self.z2Scale, self.xzScale, self.yzScale, self.x2_y2Scale, self.xyScale, self.z3Scale, self.xz2Scale, self.yz2Scale, self.zx2_zy2Scale, self.xyzScale, self.x3Scale, self.y3Scale, self.z4Scale, self.z5Scale]
        for i in range(len(lineedits)):
            lineedits[i].setText("0.0")
            scale[i].blockSignals(True)
            scale[i].setValue(0)
            scale[i].blockSignals(False)
        self.sim()

    def setValues(self, fidSurface, fwhm):
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
        self.game4Button = QtWidgets.QPushButton('All', self)
        self.game4Button.clicked.connect(lambda: self.father.startGame(4))
        grid.addWidget(self.game4Button, 4, 8)
        self.gameZButton = QtWidgets.QPushButton('Z* only', self)
        self.gameZButton.clicked.connect(lambda: self.father.startGame(zonly=True))
        grid.addWidget(self.gameZButton, 5, 8)
        self.gameStopButton = QtWidgets.QPushButton('Stop', self)
        self.gameStopButton.clicked.connect(self.father.resetGame)
        grid.addWidget(self.gameStopButton, 6, 8)
        
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
            self.father.shimFrame.setTitleText(-self.shimSim.z1Game, -self.shimSim.x1Game, -self.shimSim.y1Game, -self.shimSim.z2Game, -self.shimSim.xzGame, -self.shimSim.yzGame, -self.shimSim.x2_y2Game, -self.shimSim.xyGame, -self.shimSim.z3Game, -self.shimSim.xz2Game, -self.shimSim.yz2Game, -self.shimSim.zx2_zy2Game, -self.shimSim.xyzGame, -self.shimSim.x3Game, -self.shimSim.y3Game, -self.shimSim.z4Game, -self.shimSim.z5Game)
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
        self.shimSim.simulate(*args, **kwargs)
        self.shimPlotFrame.setData(self.shimSim.freq, self.shimSim.spectrum)
        self.shimPlotFrame.showFid()
        self.shimFrame.setValues(self.shimSim.fidSurface, self.shimSim.fwhm)

    def startGame(self, *args, **kwargs):
        self.shimSim.startGame(*args, **kwargs)
        self.shimFrame.sim()

    def resetGame(self):
        self.shimSim.resetGame()
        self.shimFrame.sim()

    def resetShims(self):
        self.shimFrame.resetShims()
        
if __name__ == '__main__':
    root = QtWidgets.QApplication(sys.argv)
    mainProgram = MainProgram(root)
    mainProgram.setWindowTitle("Shimpanzee")
    mainProgram.show()
    sys.exit(root.exec_())
