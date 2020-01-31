#!/usr/bin/env python3
import json
import pprint
from pyspdcalc import *

range = PlotRange2D(x_range=(0, 10), y_range=(20, 30), steps=(11, 11))

pp = pprint.PrettyPrinter(indent=2)
pp.pprint(range.x_range)

pp.pprint(Crystal.get_all_meta())

crystal = Crystal('KTP')

meta = crystal.get_meta()
pp.pprint(meta)

pp.pprint({ "indices": crystal.get_indices(775e-9, 273) })
