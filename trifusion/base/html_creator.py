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

title_index = 0
subcat_index = 1
img_index = 2
desc_index = 3

class HtmlTemplate:
    """
    Data must be a list of tuples (title, image, description)
    """
    def __init__(self, folder, title, data):
        if not folder.endswith("/"):
            folder += "/"
        self.folder = folder
        self.title = title
        self.data = data

    def write_file(self):
        if not self.data:
            return None

        #add header, styles basic html info
        html = "<!DOCTYPE html>\n\
<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\" />\n\
<html>\n\
<head>\n\
     <style type=\"text/css\">\n\
        .flex-container {\n\
            display: -webkit-flex;\n\
            display: -ms-flex;\n\
            display: -moz-flex;\n\
            display: flex;\n\
            -webkit-flex-direction: row;\n\
            flex-direction: row;\n\
        }\n\
        .flex-item {\n\
            display: -webkit-flex;\n\
            display: -ms-flex;\n\
            display: -moz-flex;\n\
            display: flex;\n\
            margin: 3px;\n\
            padding: 0 0 10px;\n\
        }\n\
        .flex-item img{\n\
            width: 100%;\n\
        }\n\
        span {\n\
            min-width: 8em;\n\
            padding-right: 1em;\n\
        }\n\
        a {\n\
            padding-left: 0.3em;\n\
        }\n\
        .img-wrapper {\n\
            display: inline-block;\n\
            overflow: hidden;\n\
            border: 1px solid gray;\n\
        }\n\
    </style>\n\
    <script type=\"text/javascript\">\n\
        var images = ["
        #add images names
        for (_, _, img, _) in self.data[:-1]:
            html += "\"" + img + "\","
        html += "\"" + self.data[-1][img_index] + "\"];\n"

        #add descriptions
        html += "\t\tvar descriptions = ["
        for (_, _, _, desc) in self.data[:-1]:
            html += "\"" + desc + "\","
        html += "\"" + self.data[-1][desc_index] + "\"];\n"

        #add titles
        html += "\t\tvar titles = ["
        for (title, _, _, _) in self.data[:-1]:
            html += "\"" + title + "\","
        html += "\"" + self.data[-1][title_index] + "\"];\n"

        #add data changing funtion
        html += "\t\tfunction changeImage(index){\n\
                var img = document.getElementById(\"img_place\");\n\
                img.src = \"Figures/\" + images[index];\n\
                document.getElementById(\"desc_place\").innerHTML = descriptions[index];\n\
                document.getElementById(\"subtitle\").textContent = titles[index];\n\
            }\n\
    </script>\n"

        #add actual html
        #place 1st image visible by default
        html += "</head>\n\
<body>\n\
    <h1 >" + self.title + "</h1>\
    <h2 id=\"subtitle\">" + self.data[0][0] + "</h2>\n\
    <div class=\"flex-container\">\n\
        <div class=\"flex-item\">\n\
            <span>\n"

        #for each image add code to change them
        last_cat = ""
        for x in xrange(len(self.data)):
            cat = self.data[x][subcat_index]
            if last_cat != cat:
                html += "\t\t\t\t<br><h3>" + cat + "</h3></br>\n"
                last_cat = cat
            html += "\t\t\t\t<a href=\"#\" onclick=\"changeImage(" + str(x) + ")\">" + self.data[x][title_index] + "</a></p>\n"

        #insert remainder html
        html += "\t\t\t</span>\n\
            <div class=\"img-wrapper\">\n\
                <img id=\"img_place\" src=\"Figures/" + self.data[0][img_index] + "\"/>\n\
           </div>\n\
        </div>\n\
    </div>\n\
    <h2>Description</h2>\n\
    <span id=\"desc_place\">" + self.data[0][desc_index] + "</span>\n\
</body>\n\
</html>"

        output_handle = open(self.folder + "index2.html", "w")
        output_handle.write(html)
        output_handle.close()

if __name__ == "__main__":
    html = HtmlTemplate("/home/fernando-work/Documents/trihtml", "Orthology exploration report", [("ay ay", "cat1", "Species_copy_number.png", "desc <b> cenas </b> "), ("ay ay", "cat1", "Species_copy_number.png", "desc <b> cenas </b> "), ("ay ay ay", "cat2", "Species_coverage.png", "desc2"), ("ay ay", "cat2", "Species_copy_number.png", "desc")])
    html.write_file()

__author__ = "Diogo N. Silva and Fernando Alves"
