from pysrc.lib.__init__ import *

FOLDER2SAVE = os.path.dirname(os.path.realpath(__file__)) + '/results/'
FOLDER2FIND = os.path.dirname(os.path.realpath(__file__)) + '/datafromlit/'
if not os.path.isdir(FOLDER2SAVE): os.mkdir(FOLDER2SAVE)