# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ChipPy.ui'
#
# Created: Thu Jan 23 10:27:57 2014
#      by: PyQt4 UI code generator 4.10.3
#
# WARNING! All changes made in this file will be lost!

import os, pwd, sys, ConfigParser

from PyQt4.QtCore import *
from PyQt4.QtGui import *

from chippy.util.run_record import RunRecord
from chippy.express import db_query
from chippy.util.util import run_command

try:
    _fromUtf8 = QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QApplication.translate(context, text, disambig)

class Ui_MainWindow(object):
    def buildDBGrid(self):
        """ Sets up the display area for ChipPy DB """
        self.gridLayout = QGridLayout()
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))

        # Title text
        self.title_label = QLabel(self.centralwidget)
        font = QFont()
        font.setPointSize(21)
        font.setBold(True)
        font.setWeight(75)
        self.title_label.setFont(font)
        self.title_label.setAlignment(Qt.AlignCenter)
        self.title_label.setObjectName(_fromUtf8("title_label"))
        self.gridLayout.addWidget(self.title_label, 0, 0, 1, 6)

        # standard font for fixed display text
        font = QFont()
        font.setPointSize(13)
        font.setBold(True)
        font.setWeight(75)

        # display DB path
        self.db_label = QLabel(self.centralwidget)
        self.db_label.setFont(font)
        self.db_label.setLayoutDirection(Qt.LeftToRight)
        self.db_label.setObjectName(_fromUtf8("db_label"))
        self.gridLayout.addWidget(self.db_label, 2, 0, 1, 1)
        self.current_db_label = QLabel(self.centralwidget)
        self.current_db_label.setObjectName(_fromUtf8("current_db_label"))
        self.gridLayout.addWidget(self.current_db_label, 2, 1, 1, 5)

        # Species
        self.species_label = QLabel(self.centralwidget)
        self.species_label.setFont(font)
        self.species_label.setObjectName(_fromUtf8("species_label"))
        self.gridLayout.addWidget(self.species_label, 4, 0, 1, 1)
        self.current_species_label = QLabel(self.centralwidget)
        self.current_species_label.setObjectName(_fromUtf8("current_species_label"))
        self.gridLayout.addWidget(self.current_species_label, 4, 1, 1, 1)

        # Ensembl release
        self.release_label = QLabel(self.centralwidget)
        self.release_label.setFont(font)
        self.release_label.setObjectName(_fromUtf8("release_label"))
        self.gridLayout.addWidget(self.release_label, 4, 2, 1, 1)
        self.current_release_label = QLabel(self.centralwidget)
        self.current_release_label.setObjectName(_fromUtf8("current_release_label"))
        self.gridLayout.addWidget(self.current_release_label, 4, 3, 1, 1)
        self.centralLayout.addLayout(self.gridLayout)

        # Chromosomes recognised in DB
        self.chroms_label = QLabel(self.centralwidget)
        self.chroms_label.setFont(font)
        self.chroms_label.setObjectName(_fromUtf8("chroms_label"))
        self.gridLayout.addWidget(self.chroms_label, 5, 0, 1, 1)
        self.current_chroms_label = QLabel(self.centralwidget)
        self.current_chroms_label.setObjectName(_fromUtf8("current_chroms_label"))
        self.gridLayout.addWidget(self.current_chroms_label, 5, 1, 1, 5)

        # Number of genes in DB
        self.genes_label = QLabel(self.centralwidget)
        self.genes_label.setFont(font)
        self.genes_label.setObjectName(_fromUtf8("genes_label"))
        self.gridLayout.addWidget(self.genes_label, 6, 0, 1, 1)
        self.number_genes_label = QLabel(self.centralwidget)
        self.number_genes_label.setObjectName(_fromUtf8("number_genes_label"))
        self.gridLayout.addWidget(self.number_genes_label, 6, 1, 1, 1)

        # Exons
        self.exons_label = QLabel(self.centralwidget)
        self.exons_label.setFont(font)
        self.exons_label.setObjectName(_fromUtf8("exons_label"))
        self.gridLayout.addWidget(self.exons_label, 6, 2, 1, 1)
        self.number_exons_label = QLabel(self.centralwidget)
        self.number_exons_label.setObjectName(_fromUtf8("number_exons_label"))
        self.gridLayout.addWidget(self.number_exons_label, 6, 3, 1, 1)

    def buildViewerRegion(self):
        """ QTableWidgets for DB and log file viewing """

        # Create a resizable vertical layout to hold the DB and log viewers
        self.viewerLayout = QVBoxLayout()
        self.viewerLayout.setObjectName(_fromUtf8("viewerLayout"))

        # Create a QTableWidget for DB content viewing inside a QGroupBox
        self.db_view = QGroupBox(title='Expression data in DB',
                flat=True)
        self.db_view.setObjectName(_fromUtf8('db_view'))
        self.db_view.setToolTip('Current contents of database')
        self.db_layout = QVBoxLayout(self.db_view)
        self.db_layout.setObjectName(_fromUtf8('db_layout'))
        self.db_view.setLayout(self.db_layout)
        self.viewerLayout.addWidget(self.db_view)

        # Create a table and add to the db_layout
        self.db_table = QTableWidget()
        self.db_table.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.db_table.setHorizontalScrollMode(QAbstractItemView.ScrollPerPixel)
        db_columns = ['Sample', 'Description', 'Type',
                   'Genes present', 'Reference files']
        model = QStandardItemModel(self.db_table)
        model.setHorizontalHeaderLabels(db_columns)
        header = QHeaderView(Qt.Horizontal, self.db_table)
        header.setModel(model)
        header.setStretchLastSection(True)
        header.setResizeMode(QHeaderView.ResizeToContents)
        header.setFixedHeight(25)
        self.db_table.setHorizontalHeader(header)
        self.db_layout.addWidget(self.db_table)
        self.db_table.setColumnCount(len(db_columns))

        # Create a horizontal spacer for adjusting the relative viewer sizes
        spacerItem = QSpacerItem(40, 20, QSizePolicy.Expanding,
                QSizePolicy.Minimum)
        self.viewerLayout.addItem(spacerItem)

        # Since a log file is tab delimited with fixed columns we can
        # display it as a table too
        self.log_view = QGroupBox(title='ChIPPy execution log',
                flat=True)
        self.log_view.setObjectName(_fromUtf8('log_view'))
        self.log_view.setToolTip('Record log of ChIPPy script execution')
        self.log_layout = QVBoxLayout(self.log_view)
        self.log_view.setLayout(self.log_layout)
        self.viewerLayout.addWidget(self.log_view)

        self.log_table = QTableWidget()
        self.log_table.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.log_table.setHorizontalScrollMode(QAbstractItemView.ScrollPerPixel)
        log_columns = ['Date & Time', 'Code section', 'Log level',
                      'Message', 'Value']
        model = QStandardItemModel(self.log_table)
        model.setHorizontalHeaderLabels(log_columns)
        header = QHeaderView(Qt.Horizontal, self.log_table)
        header.setModel(model)
        header.setResizeMode(QHeaderView.ResizeToContents)
        #header.setStretchLastSection(True)
        header.setFixedHeight(25)
        #header.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.log_table.setHorizontalHeader(header)
        self.log_layout.addWidget(self.log_table)
        self.log_table.setColumnCount(len(log_columns))

        # Add this block to the main layout
        self.centralLayout.addLayout(self.viewerLayout)

    def buildFileMenuMembers(self):
        """ Add a file menu and actions """
        self.menuFile = QMenu(self.menubar)
        self.menuFile.setObjectName(_fromUtf8("menuFile"))

        self.actionCreate_database = QAction(MainWindow)
        self.actionCreate_database.setObjectName(_fromUtf8("actionCreate_database"))
        self.actionCreate_database.triggered.connect(self.start_chippy_db)
        self.menuFile.addAction(self.actionCreate_database)


        self.actionOpen_database = QAction(MainWindow)
        self.actionOpen_database.setObjectName(_fromUtf8("actionOpen_database"))
        self.actionOpen_database.triggered.connect(self.open_chippy_db)
        self.menuFile.addAction(self.actionOpen_database)

        self.actionClose_database = QAction(MainWindow)
        self.actionClose_database.setObjectName(_fromUtf8("actionClose_database"))
        self.actionClose_database.triggered.connect(self.close_chippy_db)
        self.menuFile.addAction(self.actionClose_database)

        self.actionExit = QAction(MainWindow)
        self.actionExit.setObjectName(_fromUtf8("actionExit"))
        self.actionExit.setShortcut('Ctrl+Q')
        self.actionExit.setStatusTip('Exit application')
        self.actionExit.triggered.connect(qApp.quit)
        self.menuFile.addAction(self.actionExit)

        self.menubar.addAction(self.menuFile.menuAction())

    def buildFunctionalAnalysisMenuMembers(self):
        """ Add a menu and actions for functional analysis """
        self.menuFunctional_Analysis = QMenu(self.menubar)
        self.menuFunctional_Analysis.setObjectName(_fromUtf8("menuFunctional_Analysis"))
        self.menubar.addAction(self.menuFunctional_Analysis.menuAction())

        self.actionAdd_expression_data = QAction(MainWindow)
        self.actionAdd_expression_data.setObjectName(_fromUtf8("actionAdd_expression_data"))
        self.actionAdd_expression_data.triggered.connect(self.add_expression_data)
        self.actionDrop_expression_data = QAction(MainWindow)
        self.actionDrop_expression_data.setObjectName(_fromUtf8("actionDrop_expression_data"))
        self.actionDrop_expression_data.triggered.connect(self.drop_expression_data)
        self.actionExport_counts_data = QAction(MainWindow)
        self.actionExport_counts_data.setObjectName(_fromUtf8("actionExport_counts_data"))
        self.actionExport_counts_data.triggered.connect(self.export_counts_data)
        self.actionPlot_counts_v_expression = QAction(MainWindow)
        self.actionPlot_counts_v_expression.setObjectName(_fromUtf8("actionPlot_counts_v_expression"))
        self.actionPlot_counts_v_expression.triggered.connect(self.plot_counts_data)

        self.menuFunctional_Analysis.addAction(self.actionAdd_expression_data)
        self.menuFunctional_Analysis.addAction(self.actionDrop_expression_data)
        self.menuFunctional_Analysis.addAction(self.actionExport_counts_data)
        self.menuFunctional_Analysis.addAction(self.actionPlot_counts_v_expression)

    def buildDataRelationShipsMenuMembers(self):
        """ Add data relationships menu and actions """
        self.menuData_Relationships = QMenu(self.menubar)
        self.menuData_Relationships.setObjectName(_fromUtf8("menuData_Relationships"))

        self.menubar.addAction(self.menuData_Relationships.menuAction())

        self.actionExpression_difference_v_absolute = QAction(MainWindow)
        self.actionExpression_difference_v_absolute.setObjectName(_fromUtf8("actionExpression_difference_v_absolute"))
        self.actionExpression_difference_v_absolute.triggered.connect(self.expression_diff_vs_abs)
        self.actionCounts_v_Expression = QAction(MainWindow)
        self.actionCounts_v_Expression.setObjectName(_fromUtf8("actionCounts_v_Expression"))
        self.actionCounts_v_Expression.triggered.connect(self.counts_vs_expression)
        self.actionCounts_distribution = QAction(MainWindow)
        self.actionCounts_distribution.setObjectName(_fromUtf8("actionCounts_distribution"))
        self.actionCounts_distribution.triggered.connect(self.counts_distribution)
        self.actionExpression_distribution = QAction(MainWindow)
        self.actionExpression_distribution.setObjectName(_fromUtf8("actionExpression_distribution"))
        self.actionExpression_distribution.triggered.connect(self.expression_distribution)

        self.menuData_Relationships.addAction(self.actionExpression_difference_v_absolute)
        self.menuData_Relationships.addAction(self.actionCounts_v_Expression)
        self.menuData_Relationships.addAction(self.actionCounts_distribution)
        self.menuData_Relationships.addAction(self.actionExpression_distribution)

    def buildLogMenuMembers(self):
        """ add a menu and actions for the ChipPy Log """
        self.menuChipPy_Log = QMenu(self.menubar)
        self.menuChipPy_Log.setObjectName(_fromUtf8("menuChipPy_Log"))

        self.menubar.addAction(self.menuChipPy_Log.menuAction())

        self.actionEmpty_log_records = QAction(MainWindow)
        self.actionEmpty_log_records.setObjectName(_fromUtf8("actionEmpty_log_records"))
        self.actionEmpty_log_records.triggered.connect(self.empty_log_records)

        self.menuChipPy_Log.addAction(self.actionEmpty_log_records)

    def buildHelpMenuMembers(self):
        """ add help menu and actions """
        self.menuHelp = QMenu(self.menubar)
        self.menuHelp.setObjectName(_fromUtf8("menuHelp"))
        self.menubar.addAction(self.menuHelp.menuAction())

        self.actionTutorial = QAction(MainWindow)
        self.actionTutorial.setObjectName(_fromUtf8("actionTutorial"))
        self.actionTutorial.triggered.connect(self.tutorial)

        self.actionAbout = QAction(MainWindow)
        self.actionAbout.setObjectName(_fromUtf8("actionAbout"))
        self.actionAbout.triggered.connect(self.about)

        self.menuHelp.addAction(self.actionTutorial)
        self.menuHelp.addAction(self.actionAbout)

    def buildMenuBar(self):
        """ All ChipPy functionality is accessed through the menu bar """
        self.menubar = QMenuBar(MainWindow)
        self.menubar.setObjectName(_fromUtf8("menubar"))

        self.buildFileMenuMembers()
        self.buildFunctionalAnalysisMenuMembers()
        self.buildDataRelationShipsMenuMembers()
        self.buildLogMenuMembers()
        self.buildHelpMenuMembers()

        MainWindow.setMenuBar(self.menubar)

    def buildStatusBar(self):
        """ Creates a status bar - not used in this app """
        self.statusbar = QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)

    def setupUi(self, MainWindow):
        """ GUI entry point """
        self.current_db = None

        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        self.centralwidget = QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralWidget"))

        self.centralLayout = QVBoxLayout(self.centralwidget)
        self.centralLayout.setSizeConstraint(QLayout.SetDefaultConstraint)
        self.centralLayout.setMargin(10)
        self.centralLayout.setObjectName(_fromUtf8("centralLayout"))

        # Grid at the top for all the text
        self.buildDBGrid()
        self.buildViewerRegion()
        self.buildMenuBar()
        self.buildStatusBar()

        self.switch_menu_actions(False)

        #self.populateDBInfo(self.current_db_path)
        #self.populateSampleTable(self.current_db_path)
        self.populateLogTable()

        MainWindow.setCentralWidget(self.centralwidget)

        self.retranslateUi(MainWindow)
        QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        """
            This will attempt to translate labels and menu items to another
            language. It'll likely be a poor translation of these terms.
        """
        MainWindow.setWindowTitle(_translate("MainWindow", "ChIPPy: Functional Genomics Investigator", None))
        self.current_db_label.setText(_translate("MainWindow", '', None))
        self.number_exons_label.setText(_translate("MainWindow", '', None))
        self.db_label.setText(_translate("MainWindow", "Current DB", None))
        self.release_label.setText(_translate("MainWindow", "Ensembl Release", None))
        self.chroms_label.setText(_translate("MainWindow", "Chromosomes", None))
        self.genes_label.setText(_translate("MainWindow", "# Genes", None))
        self.species_label.setText(_translate("MainWindow", "Species", None))
        self.exons_label.setText(_translate("MainWindow", "# Exons", None))
        self.current_chroms_label.setText(_translate("MainWindow", '', None))
        self.number_genes_label.setText(_translate("MainWindow", '', None))
        self.title_label.setText(_translate("MainWindow", "ChIPPy: Functional Genomics Investigator", None))
        self.current_species_label.setText(_translate("MainWindow", '', None))
        self.current_release_label.setText(_translate("MainWindow", '', None))

        # Menu bar headers
        self.menuFile.setTitle(_translate("MainWindow", "File", None))
        self.menuFunctional_Analysis.setTitle(_translate("MainWindow", "Functional Analysis", None))
        self.menuData_Relationships.setTitle(_translate("MainWindow", "Data Relationships", None))
        self.menuChipPy_Log.setTitle(_translate("MainWindow", "ChipPy Log", None))
        self.menuHelp.setTitle(_translate("MainWindow", "Help", None))

        # File menu
        self.actionCreate_database.setText(_translate("MainWindow", "Create database", None))
        self.actionOpen_database.setText(_translate("MainWindow", "Open database", None))
        self.actionClose_database.setText(_translate("MainWindow", "Close database", None))
        self.actionExit.setText(_translate("MainWindow", "&Exit ChIPPy", None))

        # Functional Analysis menu
        self.actionAdd_expression_data.setText(_translate("MainWindow", "Add expression data", None))
        self.actionDrop_expression_data.setText(_translate("MainWindow", "Drop expression data", None))
        self.actionExport_counts_data.setText(_translate("MainWindow", "Export counts data", None))
        self.actionPlot_counts_v_expression.setText(_translate("MainWindow", "Plot counts v expression", None))

        # Relationships menu
        self.actionExpression_difference_v_absolute.setText(_translate("MainWindow", "Expression difference v absolute", None))
        self.actionCounts_v_Expression.setText(_translate("MainWindow", "Counts v Expression", None))
        self.actionCounts_distribution.setText(_translate("MainWindow", "Counts distribution", None))
        self.actionExpression_distribution.setText(_translate("MainWindow", "Expression distribution", None))

        # Log menu
        self.actionEmpty_log_records.setText(_translate("MainWindow", "Empty log records", None))

        # Help menu
        self.actionTutorial.setText(_translate("MainWindow", "Tutorial", None))
        self.actionAbout.setText(_translate("MainWindow", "About ChIPPy", None))

    def populateDBInfo(self, session):
        """ Display basic DB info in top panel """
        self.current_db_label.setText(self.current_db)
        chroms = sorted(db_query.get_chroms(session))
        self.current_chroms_label.setText\
            (' '.join(chroms))
        db_name_parts = self.current_db.split('_')
        release = db_name_parts[-2]
        species = db_name_parts[-1].split('.')[0]
        self.current_release_label.setText(release)
        self.current_species_label.setText(species)
        self.number_genes_label.setText(str(db_query.get_gene_counts(session)))
        self.number_exons_label.setText(str(db_query.get_exon_counts(session)))

    def populateDBTable(self, session=None):
        """ Get all expression set data from self.current_db """
        if session is None:
            if not self.check_valid_db(self.current_db):
                return
            session = db_query.make_session(self.current_db)

        names_descriptions = db_query.get_sample_names_descriptions(session)
        names = names_descriptions.keys()
        descriptions = names_descriptions.values()
        types = []
        num_genes = []
        files = []
        for name in names:
            abs = set(db_query.get_expr_sample_names(session))
            diff = set(db_query.get_diff_sample_names(session))
            target = set(db_query.get_target_gene_names(session))

            if name in abs:
                types.append('Expression')
                num_genes.append(db_query.get_expression_counts(session,
                        sample_name=name))
            elif name in diff:
                types.append('Differential')
                num_genes.append(db_query.get_diff_counts(session,
                        sample_name=name))
            elif name in target:
                types.append('Target Genes')
                num_genes.append(db_query.get_targetgene_counts(session,
                        sample_name=name))
            # reffile_entries returns reffile objects
            reffiles = db_query.get_reffile_entries(session,
                sample_name=name)
            file_names = [r.name for r in reffiles]
            files.append(', '.join(file_names))
        session.close()

        self.db_table.setRowCount(0)
        for row, (n, d, t, g, f) in enumerate(zip(names, descriptions, types, num_genes, files)):
            self.db_table.setRowCount(self.db_table.rowCount()+1)
            self.db_table.setItem(row, 0, QTableWidgetItem(QString(n)))
            self.db_table.setItem(row, 1, QTableWidgetItem(QString(d)))
            self.db_table.setItem(row, 2, QTableWidgetItem(QString(t)))
            self.db_table.setItem(row, 3, QTableWidgetItem(QString(str(g))))
            self.db_table.setItem(row, 4, QTableWidgetItem(QString(f)))

    def populateLogTable(self):
        """ display ChipPy Log text in the appropriate window """
        rr = RunRecord()
        self.log_table.setRowCount(0)
        try:
            table = rr.getMessageTable(last_n_lines=30, include_date=True)
        except RuntimeError:
            return

        if table is None:
            return
        else:
            for r, row in enumerate(table):
                self.log_table.setRowCount(self.log_table.rowCount()+1)
                for c, column in enumerate(row):
                    self.log_table.setItem(r,c,QTableWidgetItem(QString(column)))

    # file menu action functions

    def switch_menu_actions(self, enable):
        """ enable on valid DB, disable if no DB """
        self.actionAdd_expression_data.setEnabled(enable)
        self.actionDrop_expression_data.setEnabled(enable)
        self.actionExport_counts_data.setEnabled(enable)
        self.actionPlot_counts_v_expression.setEnabled(enable)

        self.actionExpression_difference_v_absolute.setEnabled(enable)
        self.actionCounts_v_Expression.setEnabled(enable)
        self.actionCounts_distribution.setEnabled(enable)
        self.actionExpression_distribution.setEnabled(enable)

    def check_valid_db(self, db_path):
        """ True if valid data in DB at path """
        if db_path is None or db_path == '':
            return False

        # test DB is valid
        session = db_query.make_session(db_path)
        if db_query.get_species(session) is None:
            session.close()
            return False

        session.close()
        return True

    def start_chippy_db(self):
        """ script to create a new DB """
        command = self._make_cmd_str('start_chippy_db.py', include_db=False)
        returncode, stdout, stderr = run_command(command)
        if returncode == 0:
            if self.check_valid_db(stdout):
                self.current_db = stdout
                # Check the DB works correctly
                session = db_query.make_session(self.current_db)
                self.populateDBInfo(session)
                session.close()
                self.switch_menu_actions(True)

    def open_chippy_db(self):
        """ Use dialog to select DB file and populate view with DB info """
        rr = RunRecord('open_chippy_db')
        db_path = str(QFileDialog.getOpenFileName())
        if not self.check_valid_db(db_path):
            rr.addWarning('DB has invalid format', db_path)
            self.populateLogTable()
            return
        self.current_db = os.path.realpath(db_path)

        session = db_query.make_session(self.current_db)
        self.populateDBInfo(session)
        self.populateDBTable(session)
        rr.addInfo('DB opened successfully', db_path)
        self.populateLogTable()
        session.close()
        self.switch_menu_actions(True)

    def close_chippy_db(self):
        """ Set DB to None and clear DB display text """
        self.current_db = None
        self.db_table.setRowCount(0)
        self.current_chroms_label.setText('')
        self.current_db_label.setText('')
        self.current_release_label.setText('')
        self.current_species_label.setText('')
        self.number_genes_label.setText('')
        self.number_exons_label.setText('')

        self.switch_menu_actions(False)

    # functional analysis menu action functions

    def _make_cmd_str(self, script_name, include_db=True):
        """ build command line for launching scripts """
        if include_db:
            return 'cd scripts; python ' + \
                   script_name + ' ' + self.current_db
        else:
            return 'cd scripts; python ' + script_name

    def add_expression_data(self):
        """ launch add_expression_db.py """
        command = self._make_cmd_str('add_expression_db.py')
        run_command(command)
        self.populateDBTable()
        self.populateLogTable()

    def drop_expression_data(self):
        """ launch drop_expression_db.py """
        command = self._make_cmd_str('drop_expression_db.py')
        run_command(command)
        self.populateDBTable()
        self.populateLogTable()

    def export_counts_data(self):
        """ launch export_counts.py """
        command = self._make_cmd_str('export_counts.py')
        run_command(command)
        self.populateLogTable()

    def plot_counts_data(self):
        """ launch plot_counts.py """
        command = self._make_cmd_str('plot_counts.py')
        run_command(command)
        self.populateLogTable()

    # data relationships menu action functions

    def expression_diff_vs_abs(self):
        """ launch diff_abs_plots.py """
        command = self._make_cmd_str('diff_abs_plots.py')
        run_command(command)
        self.populateLogTable()

    def counts_vs_expression(self):
        """ launch counts_vs_expr.py """
        command = self._make_cmd_str('counts_vs_expr.py')
        run_command(command)
        self.populateLogTable()

    def counts_distribution(self):
        """ launch counts_distribution.py """
        command = self._make_cmd_str('counts_distribution.py')
        run_command(command)
        self.populateLogTable()

    def expression_distribution(self):
        """ launch expr_distribution.py """
        command = self._make_cmd_str('expr_distribution.py')
        run_command(command)
        self.populateLogTable()

    # log menu action functions

    def empty_log_records(self):
        """ call run_record.clear_log() """
        rr = RunRecord()
        rr.empty_log()
        self.populateLogTable()

    def tutorial(self):
        """ Display a usage tutorial to the user """
        #TODO: write a tutorial and display it!
        #TODO: convert to HTML
        tut_dialog = QDialog(self.centralwidget)
        layout = QVBoxLayout(tut_dialog)
        text_edit = QTextBrowser()
        with open('ChipPy_manual.txt') as f:
            contents = f.read()
        text_edit.setText(contents)
        layout.addWidget(text_edit)

        tut_dialog.setLayout(layout)
        tut_dialog.setGeometry(300, 300, 350, 300)
        tut_dialog.setWindowTitle('ChipPy Manual and Tutorial')
        tut_dialog.show()
        tut_dialog.raise_()
        tut_dialog.activateWindow()

    def about(self):
        """ Display authors and version info """
        title_txt = QString('ChipPy: Functional Genomics Investigator')
        authors_txt = QString('Cameron A. Jack, Anuj Pahwa, Gavin A. '+\
            'Huttley, 2010-2014. The John Curtin School of Medical '+\
            'Research, The Australian National University.')
        info_txt = QString('ChipPy shows the relationship between '+\
            'mapped targets and gene expression at annotated gene feature '+\
            'locations.')
        dialog = QMessageBox()
        QMessageBox.information(dialog, title_txt, authors_txt+'<br>'+info_txt)

if __name__ == "__main__":
    config = ConfigParser.ConfigParser()
    with open('scripts/chippy.ini', 'w') as cfg_file:
        config.add_section('Log')
        config.set('Log', 'Directory', os.path.realpath(os.curdir))
        config.write(cfg_file)
        cfg_file.close()

    rr = RunRecord()
    rr.addInfo('ChipPy launched by', pwd.getpwuid(os.getuid()).pw_name)
    app = QApplication(sys.argv)
    MainWindow = QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

