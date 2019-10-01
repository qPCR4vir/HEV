from pathlib import Path
import logging

import sys
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QMainWindow, QAction,     qApp,        QApplication, QFileDialog, \
                            QTableView,  QVBoxLayout, QMessageBox
from PyQt5.QtGui import QIcon
from PyQt5 import QtSql

print('import seq_sql')
from seq_sql import create_all


class MW(QMainWindow):
    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):

        self.setGeometry(300, 100, 1300, 800)
        self.setWindowTitle('Genotyping')

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

        origData = QAction('View &Origial Data', self)
        origData.triggered.connect(self.viewOrigData)

        self.statusBar()

        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(newProjAct)
        fileMenu.addAction(openProjAct)
        fileMenu.addAction(closeProjAct)
        fileMenu.addAction(GBfileAct)
        fileMenu.addAction(impAlignAct)
        fileMenu.addAction(exitAct)

        viewMenu = menubar.addMenu('&View')
        # viewMenu.
        viewMenu.addAction(origData)

        self.toolbar = self.addToolBar('main')
        self.toolbar.addAction(newProjAct)
        self.toolbar.addAction(openProjAct)
        self.toolbar.addAction(closeProjAct)
        self.toolbar.addAction(GBfileAct)
        self.toolbar.addAction(impAlignAct)
        self.toolbar.addAction(exitAct)

        self.v = QTableView()
        self.setCentralWidget(self.v)

        #vbox = QVBoxLayout()
        #vbox.addWidget(self.v)

        #self.setLayout(vbox)


        self.show()

    def newDB(self):
        create_all()
        return

    def openDB(self):
        print('opening...')
        dbfilename = QFileDialog.getOpenFileName(self, 'Select the sequence databank to open', '../data/temp/seq- (30).db')
        print(dbfilename)
        if not dbfilename:
            QMessageBox.critical(None, "Cannot open database",  "No fileneme specified.\n.\n\n"
                                 "Click Cancel to return.",       QMessageBox.Cancel)

        db = QtSql.QSqlDatabase.addDatabase('QSQLITE')
        db.setDatabaseName(dbfilename[0])

        if not db.open():
            QMessageBox.critical(None, "Cannot open database",
                                 "Unable to establish a database connection.\n"
                                 "This example needs SQLite support. Please read the Qt SQL "
                                 "driver documentation for information how to build it.\n\n"
                                 "Click Cancel to exit.",
                                 QMessageBox.Cancel)
            return False
        #print('SELECT * FROM original_data')
        #query = QtSql.QSqlQuery()
        #query.prepare('SELECT * FROM original_data')
        #for r in query.exec_( ):
        #    print(r)

        self.table = QtSql.QSqlTableModel(db=db)
        self.table.setTable('original_data')
        self.table.setEditStrategy(QtSql.QSqlTableModel.OnManualSubmit)
        self.table.select()

        self.table.setHeaderData(0, Qt.Horizontal, 'Strain')
        self.table.setHeaderData(1, Qt.Horizontal, 'genotype')
        self.table.setHeaderData(2, Qt.Horizontal, "Isolate")
        #rc = self.table.rowCount()
        #print('Row count =', rc)
        #self.table.selectRow(self.table.rowCount()-1)

        self.v.setModel(self.table)
        #self.v.
        self.v.resizeColumnsToContents()
        f = self.v.verticalHeader().font()
        f.setPixelSize(8)
        #self.v.selectRow(self.table.rowCount()-1)
        #m = self.v.verticalHeader().getContentsMargins()
        #print('getContentsMargins' , m)
        #self.v.verticalHeader().
        self.v.setSortingEnabled(True)
        while self.table.canFetchMore():
            self.table.fetchMore()
        self.v.verticalHeader().setFont(f)
        self.v.verticalHeader().setDefaultSectionSize(12)
        #self.v.verticalHeader().setDefaultSectionSize(self.v.verticalHeader().minimumSectionSize())

        return

    def viewOrigData(self):
        self.table.setTable('original_data')
        self.table.setEditStrategy(QtSql.QSqlTableModel.OnManualSubmit)
        self.table.select()

        self.table.setHeaderData(0, Qt.Horizontal, 'Strain')
        self.table.setHeaderData(1, Qt.Horizontal, 'genotype')
        self.table.setHeaderData(2, Qt.Horizontal, "Isolate")
        #rc = self.table.rowCount()
        #print('Row count =', rc)
        #self.table.selectRow(self.table.rowCount()-1)

        self.v.setModel(self.table)
        #self.v.
        self.v.resizeColumnsToContents()
        f = self.v.verticalHeader().font()
        f.setPixelSize(8)
        #self.v.selectRow(self.table.rowCount()-1)
        #m = self.v.verticalHeader().getContentsMargins()
        #print('getContentsMargins' , m)
        #self.v.verticalHeader().
        self.v.setSortingEnabled(True)
        while self.table.canFetchMore():
            self.table.fetchMore()
        self.v.verticalHeader().setFont(f)
        self.v.verticalHeader().setDefaultSectionSize(12)
        #self.v.verticalHeader().setDefaultSectionSize(self.v.verticalHeader().minimumSectionSize())
        return

    def viewStrains(self):
        self.table.setTable('original_data')
        self.table.setEditStrategy(QtSql.QSqlTableModel.OnManualSubmit)
        self.table.select()

        self.table.setHeaderData(0, Qt.Horizontal, 'Strain')
        self.table.setHeaderData(1, Qt.Horizontal, 'genotype')
        self.table.setHeaderData(2, Qt.Horizontal, "Isolate")
        #rc = self.table.rowCount()
        #print('Row count =', rc)
        #self.table.selectRow(self.table.rowCount()-1)

        self.v.setModel(self.table)
        #self.v.
        self.v.resizeColumnsToContents()
        f = self.v.verticalHeader().font()
        f.setPixelSize(8)
        #self.v.selectRow(self.table.rowCount()-1)
        #m = self.v.verticalHeader().getContentsMargins()
        #print('getContentsMargins' , m)
        #self.v.verticalHeader().
        self.v.setSortingEnabled(True)
        while self.table.canFetchMore():
            self.table.fetchMore()
        self.v.verticalHeader().setFont(f)
        self.v.verticalHeader().setDefaultSectionSize(12)
        #self.v.verticalHeader().setDefaultSectionSize(self.v.verticalHeader().minimumSectionSize())
        return


    def closeDB(self):
        return

    def importGB(self):
        return


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = MW()
    sys.exit(app.exec_())
