#!/usr/bin/env python2
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


class HtmlTemplate ():
    def __init__(self):
        self.hcontents = []
        self.bcontents = []

    def add_title(self, title):
        htitle = "<title> %s </title>\n" % title
        self.hcontents.extend([htitle])

    def add_text(self, text):
        """ Adds plain text to body (no tags) """
        content = text
        self.bcontents.extend([content])

    def add_single_plot(self, heading, plot_file, heading_level="1"):
        """ Adds single plot with heading """
        head = "<h%s> %s </h1>\n" % (heading_level, heading)
        plot = "<figure> <embed type='image/svg+xml' src='%s' /> </figure>\n" %\
               plot_file
        self.bcontents.extend([head, plot])

    #def add_plot_grid(self):
        #""" TODO """

    def write_file(self, file_name):
        string = """<!DOCTYPE html>
<html>
    <head>
    %s
    </head>
    <body style='background-color:#eeeeee'>
    <h1> Alignments Report </h1>
    <p> Alignment report results ran on #DATE using options #OPTIONS
    %s
    </body>
</html>

        """ % ("\n".join(self.hcontents), "\n".join(self.bcontents))

        output_handle = open(file_name + ".html", "w")
        output_handle.write(string)
        output_handle.close()
