#!/usr/bin/python2
# -*- coding: utf-8 -*-

import sys

reload(sys)
sys.setdefaultencoding("utf8")

def gui_exec():

    try:
        from trifusion.app import main
    except ImportError:
        from app import main

    main()

if __name__ == "__main__":

    import multiprocessing

    multiprocessing.freeze_support()

    gui_exec()
