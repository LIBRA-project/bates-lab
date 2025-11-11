import openmc
from openmc.model import RectangularParallelepiped as RPP
import numpy as np
import os
from pathlib import Path
from bates_bare_model import (build_bates_model, 
                              BPE, lead, PortlandConc, steel,
                              Air,
                              concrete_wall_rpp)

source_room = 'warehouse'
source_side = 'south'
# version_string = 'sides_24in_conc_east_12in_BPE_ceil_12in_BPE_v1'
version_string = 'sides_24in_conc_no_east_BPE_ceil_18in_BPE_v1'
num_particles_per_batch = 5e5
directory = Path(source_room) / version_string / source_side

concrete_wall_mid_y = (concrete_wall_rpp.ymin.y0 + concrete_wall_rpp.ymax.y0)/2

west_shield_rpp = RPP(concrete_wall_rpp.xmin.x0 - 400,
                      concrete_wall_rpp.xmin.x0 - 400 + 24 * 2.54,
                      concrete_wall_mid_y - 200,
                      concrete_wall_mid_y + 200,
                      concrete_wall_rpp.zmin.z0,
                      concrete_wall_rpp.zmin.z0 + 8.0 * 2.54 * 12)
north_shield_rpp = RPP(west_shield_rpp.xmax.x0,
                       west_shield_rpp.xmax.x0 + 6 * 12 * 2.54,
                       west_shield_rpp.ymax.y0 - 24 * 2.54,
                       west_shield_rpp.ymax.y0,
                       west_shield_rpp.zmin.z0,
                       west_shield_rpp.zmax.z0)
south_shield_rpp = RPP(west_shield_rpp.xmax.x0,
                       west_shield_rpp.xmax.x0 + 6 * 12 * 2.54,
                       west_shield_rpp.ymin.y0,
                       west_shield_rpp.ymin.y0 + 24 * 2.54,
                       west_shield_rpp.zmin.z0,
                       west_shield_rpp.zmax.z0)
middle_shield_rpp = RPP(west_shield_rpp.xmax.x0,
                        west_shield_rpp.xmax.x0 + 6 * 12 * 2.54,
                       (north_shield_rpp.ymin.y0 + south_shield_rpp.ymax.y0)/2 - 6 * 2.54,
                       (north_shield_rpp.ymin.y0 + south_shield_rpp.ymax.y0)/2 + 6 * 2.54,
                       west_shield_rpp.zmin.z0,
                       west_shield_rpp.zmax.z0)

east_shield_rpp = RPP(south_shield_rpp.xmax.x0,
                      south_shield_rpp.xmax.x0 + 0.1 * 2.54,
                      south_shield_rpp.ymin.y0,
                      north_shield_rpp.ymax.y0,
                      west_shield_rpp.zmin.z0,
                      west_shield_rpp.zmax.z0)
ceiling_plate_rpp = RPP(west_shield_rpp.xmin.x0,
                        east_shield_rpp.xmax.x0,
                        south_shield_rpp.ymin.y0,
                        north_shield_rpp.ymax.y0,
                        west_shield_rpp.zmax.z0,
                        west_shield_rpp.zmax.z0 + 1 * 2.54)
ceiling_shield_rpp = RPP(ceiling_plate_rpp.xmin.x0,
                         ceiling_plate_rpp.xmax.x0,
                         ceiling_plate_rpp.ymin.y0,
                         ceiling_plate_rpp.ymax.y0,
                         ceiling_plate_rpp.zmax.z0,
                         ceiling_plate_rpp.zmax.z0 + 18 * 2.54)

bounding_rpp = RPP(concrete_wall_rpp.xmin.x0 - 1000,
                   concrete_wall_rpp.xmin.x0 + 100,
                   concrete_wall_rpp.ymin.y0 - 400,
                   concrete_wall_rpp.ymax.y0 + 400,
                   0,
                   400)

if source_side == 'north':
    source_point = [(west_shield_rpp.xmax.x0 + east_shield_rpp.xmin.x0)/2,
                    (north_shield_rpp.ymin.y0 + middle_shield_rpp.ymax.y0)/2,
                    50]
elif source_side == 'south':
    source_point = [(west_shield_rpp.xmax.x0 + east_shield_rpp.xmin.x0)/2,
                    (south_shield_rpp.ymax.y0 + middle_shield_rpp.ymin.y0)/2,
                    50]
else:
    raise ValueError('Invalid source side')


# west_shield_rpp = RPP(-4.0 * 2.54 * 12,
#                     -2.0 * 2.54 * 12,
#                     -4.0 * 2.54 * 12,
#                     4.0 * 2.54 * 12,
#                     0.0,
#                     8.0 * 2.54 * 12)
# north_shield_rpp = RPP(-2.0 * 2.54 * 12,
#                        3.0 * 2.54 * 12,
#                        2.0 * 2.54 * 12,
#                        4.0 * 2.54 * 12,
#                        0.0,
#                        8.0 * 2.54 * 12)
# south_shield_rpp = RPP(-2.0 * 2.54 * 12,
#                        3.0 * 2.54 * 12,
#                        -4.0 * 2.54 * 12,
#                        -2.0 * 2.54 * 12,
#                        0.0,
#                        8.0 * 2.54 * 12)
# east_shield_1_rpp = RPP(6.0 * 2.54 * 12,
#                         8.0 * 2.54 * 12,
#                         -4.0 * 2.54 * 12,
#                         4.0 * 2.54 * 12,
#                         0.0,
#                         8.0 * 2.54 * 12)

# ceiling_plate_rpp = RPP(-4.0 * 2.54 * 12,
#                         3.0 * 2.54 * 12,
#                         -4.0 * 2.54 * 12,
#                         4.0 * 2.54 * 12,
#                         96 * 2.54,
#                         97 * 2.54)


# ceiling_shield_rpp = RPP(-4.0 * 2.54 * 12,
#                          3.0 * 2.54 * 12,
#                         -4.0 * 2.54 * 12,
#                         4.0 * 2.54 * 12,
#                         97 * 2.54,
#                         109 * 2.54)

# bounding_rpp = RPP(-500, 500, - 400, 400, 0, 400)

experiment_air_region = -bounding_rpp & +west_shield_rpp & +north_shield_rpp \
    & +south_shield_rpp & +east_shield_rpp & +ceiling_shield_rpp & +ceiling_plate_rpp & +middle_shield_rpp

west_shield_cell = openmc.Cell(region=-west_shield_rpp, fill=PortlandConc, name='West Shield')
north_shield_cell = openmc.Cell(region=-north_shield_rpp, fill=PortlandConc, name='North Shield')
middle_shield_cell = openmc.Cell(region=-middle_shield_rpp, fill=BPE, name='Middle Shield')
south_shield_cell = openmc.Cell(region=-south_shield_rpp, fill=PortlandConc, name='South Shield')
east_shield_cell = openmc.Cell(region=-east_shield_rpp, fill=Air, name='East Shield 1')
ceiling_shield_cell = openmc.Cell(region=-ceiling_shield_rpp, fill=BPE, name='Ceiling Shield')
ceiling_plate_cell = openmc.Cell(region=-ceiling_plate_rpp, fill=steel, name='Ceiling Plate')
experiment_air_cell = openmc.Cell(region=experiment_air_region, fill=Air, name='Experiment Air')

universe = openmc.Universe(cells=[west_shield_cell,
                               north_shield_cell,
                               middle_shield_cell,
                               south_shield_cell,
                               east_shield_cell,
                               ceiling_shield_cell,
                               ceiling_plate_cell,
                               experiment_air_cell])

print("Source Point: ", source_point)
model = build_bates_model(experiment_universe=universe,
                          source_room=source_room,
                          num_particles_per_batch=num_particles_per_batch,
                          do_translation=False,
                          libra_center_coord=source_point)

if __name__ == '__main__':
    curr_dir = os.getcwd()
    os.makedirs(directory, exist_ok=True)
    os.chdir(directory)
    model.export_to_model_xml()
    model.plot_geometry()
    model.run(threads=14)
    os.chdir(curr_dir)


