#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#
#  Copyright 2012 Unknown <diogo@arch>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#  Author: Diogo N. Silva
#  Version: 0.1
#  Last update: 11/02/14

import mysql.connector
from os.path import join

################################################################################

"""
db_tools handles the connection and modification of MySQL databases for the
orthoMCL program
"""

sql_user = "orthouser_id"
sql_user_pass = "orthouser"
mysql_db = "orthomcl_db"


def create_orthomcl_cfg(work_dir):

    cfg_text_template = """dbVendor=mysql
dbConnectString=dbi:mysql:%s
dbLogin=%s
dbPassword=%s
similarSequencesTable=SimilarSequences
orthologTable=Ortholog1
inParalogTable=InParalog1
coOrthologTable=CoOrtholog1
interTaxonMatchView=InterTaxonMatch
percentMatchCutoff=50
evalueExponentCutoff=-5
oracleIndexTblSpc=NONE
""" % (mysql_db, sql_user, sql_user_pass)

    cfg_handle = open(join(work_dir, "orthomcl.config"), "w")
    cfg_handle.write(cfg_text_template)
    cfg_handle.close()


def sql_setup(sql_pass):

    # Connect to mysql database and handles any exception
    try:
        cnx = mysql.connector.connect(user="root", password=sql_pass,
                                     host="localhost")

        cursor = cnx.cursor()

        # Create new OrthoMCL table, drop old version if exists
        cursor.execute("DROP DATABASE IF EXISTS %s" % mysql_db)
        cursor.execute("CREATE DATABASE %s" % mysql_db)

        # Create user even if it does not exist
        cursor.execute("GRANT USAGE ON *.* TO '%s'@'localhost' IDENTIFIED BY "
                       "'%s'" % (sql_user, sql_user_pass))

        # Grant permissions to user mysql database
        cursor.execute("GRANT ALL PRIVILEGES ON %s.* TO %s@localhost IDENTIFIED"
                       " BY '%s'" % (mysql_db, sql_user, sql_user_pass))

        # Close connection
        cursor.close()
        cnx.close()

    except mysql.connector.Error as err:
        return err

__author__ = 'diogo'