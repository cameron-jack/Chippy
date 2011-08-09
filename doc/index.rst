############
Using ChipPy
############

Assuming you've been given the directory (or directories) containing the results of mapping the sequence reads to a genome, then the basic steps for using ChipPy to produce plots of tag count distributions are described below.

************************************************
Create the DB for a specific release of Ensembl.
************************************************

If you haven't been already given this file, or downloaded it our assembla site, then you need to create a new ``chippy.db`` file.

********************************************
Set your ``CHIPPY_DB`` environment variable.
********************************************

For PyCogent.app, you add this variable to the Python preferences. Specifically, open the app preferences and set a System path variable to point to where you store the ``chippy.db`` file.

For command line use, define and export this environment variable in your ``~/.bashrc``.

******************************************
Add expression data sets into the database
******************************************

Choose the file containing the results of expression analysis (exported from R in our case) for import. The data are added to the ``chippy.db`` file. You can add multiple such files to the database. These only need to be added once.

*************************
Export TSS centred counts
*************************

Choose the expression data you want to order genes by and the directory containing the mapped reads. This creates a compressed archive of counts that you use for plotting.

***********************
Plot the centred counts
***********************

Takes the results of the previous step and plots them. You can either specify the filename to save to, or simply have it draw to screen.

.. note:: Drawing hundreds of lines with thousands of data points to screen can take a while!

*************************************
PyCogent.app Bundle/ChipPy menu items
*************************************

.. csv-table::
    :header: Menu item, Task
    :widths: 5, 10

    Initialise Db, Starts a new database (Db)
    Add expression data, Adds ranked expression data from a microarray analysis to the Db
    Delete expression data, Use this to delete expression data from the Db
    Db report, Lists the files the expression data derived from
    Export centred counts, Saves tag counts centred on the TSS for expressed genes
    Plot centred counts, Takes results of above and plots them

