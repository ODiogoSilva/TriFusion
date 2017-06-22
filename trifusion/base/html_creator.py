#!/usr/bin/env python2
# -*- coding: utf-8 -*-
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


class HtmlTemplate(object):
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
            width: 20em;\n\
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
        .body{\n\
            background-color: #A0A0A0\n\
        }\n\
        h1{\n\
            color: white;\n\
            padding-top: 1em;\n\
        }\n\
        .bootstrap {\n\
            width: fill-parent;\n\
            background-color: #000000;\n\
        }\n\
        .dropdown {\n\
            position: relative;\n\
            display: inline-block;\n\
            color: white;\n\
            padding-right: 2em;\n\
            padding-left: 0.5em;\n\
        }\n\
        .dropdown-content {\n\
            display: none;\n\
            position: absolute;\n\
            background-color: #FF0000;\n\
            min-width: 15em;\n\
            max-width: 30em;\n\
            padding: 12px 16px;\n\
        }\n\
        .dropdown:hover .dropdown-content {\n\
            display: block;\n\
            color: black;\n\
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
        window.zoomedIn = false;\n\
        function cenas() {\n\
            var el = this, elp = document.getElementById(\"zoom-container\");\n\
            var zoomContainer = document.getElementById(\"img_wrapper\");\n\
            if (window.zoomedIn) {\n\
                elp.setAttribute(\"style\", \"overflow: auto\");\n\
                zoomContainer.setAttribute(\"style\", \"transform :\");\n\
                window.zoomedIn = false;\n\
            } else {\n\
                var top = el.offsetTop;\n\
                var left = el.offsetLeft - 0.25*zoomContainer.clientWidth;\n\
                var tro = (Math.abs(elp.offsetTop - el.offsetTop) > 0) ? \"bottom\" : \"top\";\n\
                tro += (Math.abs(elp.offsetLeft - el.offsetLeft) > 0) ? \" right\" : \" left\";\n\
                zoomContainer.setAttribute(\"style\", \"transform-origin: \"+ tro + \" 0px; transform: scale(2);\");\n\
                window.zoomedIn = true;\n\
            }\n\
        }\n\
    </script>\n"

        #add actual html
        #place 1st image visible by default
        html += "</head>\n\
<body>\n\
    <div class=\"bootstrap\">\n\
        <h1>" + self.title + "</h1>\n"

        #for each image add code to change them
        last_cat = None
        for x in xrange(len(self.data)):
            cat = self.data[x][subcat_index]
            if last_cat != cat:
                if last_cat:
                    html += "\t\t\t</div>\n"
                    html += "\t\t</div>\n"
                html += "\t\t<div class=\"dropdown\">\n"
                html += "\t\t\t<h3>" + cat + "</h3>\n"
                html += "\t\t\t<div class=\"dropdown-content\">\n"
                last_cat = cat
            html += "\t\t\t\t<a href=\"#\" onclick=\"changeImage(" + str(x) + ")\">" + self.data[x][title_index] + "</a></p>\n"

        html += "\t\t\t</div>\n"
        html += "\t\t</div>\n"
        html += "\t</div>\n"
        html += "\t<h2 id=\"subtitle\">" + self.data[0][0] + "</h2>\n\
    <div class=\"flex-container\">\n\
        <div class=\"flex-item\">\n\
            <div id=\"img_wrapper\" class=\"img-wrapper\">\n\
                <img id=\"img_place\" src=\"Figures/" + self.data[0][img_index] + "\" onclick=\"cenas()\" />\n\
            </div>\n\
        </div>\n\
    </div>\n\
    <h2>Description</h2>\n\
    <span id=\"desc_place\">" + self.data[0][desc_index] + "\
    </span>\n\
</body>\n\
</html>"

        output_handle = open(self.folder + "index.html", "w")
        output_handle.write(html)
        output_handle.close()

if __name__ == "__main__":
    html = HtmlTemplate("/home/fernando/Documents/trihtml", "Orthology exploration report", [("Distribution of maximum gene copies", "cat1", "Species_copy_number.png", " Derived allele frequency"), ("Species distribution", "cat2", "Species_coverage.png", "desc2"), ("Data coverage per species", "cat2", "Species_copy_number.png", "desc")])
    html.write_file()

__author__ = "Diogo N. Silva and Fernando Alves"
