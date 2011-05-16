############
Using ChipPy
############

Note that many of the menu items require you have ``CHIPPY_DB`` environment variable already set. You do this through PyCogent.app Python preferences. Specifically, by setting a System path variable to point to where you store the ``chippy.db`` file. (This file is either given to you or you've created it by running the Initialise Db step.)

.. csv-table:: Available script forms
    :header: Menu item, Task
    :widths: 5, 10
    
    Initialise Db, Starts a new database (Db)
    Add expression data, Adds ranked expression data from a microarray analysis to the Db
    Delete expression data, Use this to delete expression data from the Db
    Db report, Lists the files the expression data derived from
    Export centred counts, Saves tag counts centred on the TSS for expressed genes
    Plot centred counts, Takes results of above and plots them

