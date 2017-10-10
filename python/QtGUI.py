
import sys
from PyQt5.QtWidgets import QMainWindow, QAction, qApp, QApplication
from PyQt5.QtGui import QIcon

print('import seq_sql')
from seq_sql import create_all




class MW(QMainWindow):
    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):
        exitAct = QAction(QIcon('exit.png'), '&Exit', self)
        exitAct.setShortcut('Ctrl+Q')
        exitAct.setStatusTip('Exit application')
        exitAct.triggered.connect(qApp.quit)

        newProjAct = QAction(QIcon('database.png'), '&New virus DB project', self)
        newProjAct.setShortcut('Ctrl+N')
        newProjAct.setStatusTip('Create a new empty data bank project for virus sequences')
        newProjAct.triggered.connect(self.newDB)

        openProjAct = QAction(QIcon('open.png'), '&Open project', self)
        openProjAct.setShortcut('Ctrl+O')
        openProjAct.setStatusTip('Open an existing data bank project')
        openProjAct.triggered.connect(self.openDB)

        closeProjAct = QAction(QIcon('close.png'), '&Close project', self)
        closeProjAct.setShortcut('Ctrl+C')
        closeProjAct.setStatusTip('Close the data bank project')
        closeProjAct.triggered.connect(self.openDB)

        GBfileAct = QAction(QIcon('GB.png'), 'Import &GenBank', self)
        GBfileAct.setShortcut('Ctrl+I')
        GBfileAct.setStatusTip('Import sequence data from a flat Genbank file')
        GBfileAct.triggered.connect(self.openDB)

        impAlignAct = QAction(QIcon('fasta.png'), '&Import alignment', self)
        impAlignAct.setShortcut('Ctrl+A')
        impAlignAct.setStatusTip('Import a sequence alignment file in fasta format')
        impAlignAct.triggered.connect(self.openDB)

        self.statusBar()

        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(newProjAct)
        fileMenu.addAction(openProjAct)
        fileMenu.addAction(closeProjAct)
        fileMenu.addAction(GBfileAct)
        fileMenu.addAction(impAlignAct)
        fileMenu.addAction(exitAct)

        self.toolbar = self.addToolBar('main')
        self.toolbar.addAction(newProjAct)
        self.toolbar.addAction(openProjAct)
        self.toolbar.addAction(closeProjAct)
        self.toolbar.addAction(GBfileAct)
        self.toolbar.addAction(impAlignAct)
        self.toolbar.addAction(exitAct)


        self.setGeometry(300, 100, 1300, 800)
        self.setWindowTitle('Genotyping')
        self.show()

    def newDB(self):
        create_all()

    def openDB(self):
        pass

    def closeDB(self):
        pass

    def importGB(self):
        pass




if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = MW()
    sys.exit(app.exec_())