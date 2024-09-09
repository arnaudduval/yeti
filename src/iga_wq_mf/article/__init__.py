from pysrc.lib.__init__ import *

FOLDER2SAVE = os.path.dirname(os.path.realpath(__file__)) + '/results/'
if not os.path.isdir(FOLDER2SAVE): os.mkdir(FOLDER2SAVE)

FOLDER2DATA = os.path.dirname(os.path.realpath(__file__)) + '/datafromsimu/'
if not os.path.isdir(FOLDER2DATA): os.mkdir(FOLDER2DATA)
