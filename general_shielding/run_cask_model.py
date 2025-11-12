import openmc
from openmc.model import RectangularParallelepiped as RPP
import numpy as np
import os
from pathlib import Path
from bates_bare_model import (build_bates_model, 
                              BPE, lead, PortlandConc, steel,
                              Air)

source_room = 'cask'
source_side = 'east'
version_string = 'sides_12in_BPE_north_16in_BPE_ceil_12in_BPE_v1'
num_particles_per_batch = 1e6
directory = Path(source_room) / version_string / source_side

experiment_rpp = RPP(826.0, 1477.0, 2076.0, 2410.0, 0.0, 240.0)


west_shield_rpp = RPP(experiment_rpp.xmin.x0 + 100, 
                      experiment_rpp.xmin.x0 + 100 + 12 * 2.54,
                      experiment_rpp.ymin.y0,
                      experiment_rpp.ymin.y0 + 4 * 2.54 * 12,
                      experiment_rpp.zmin.z0,
                      experiment_rpp.zmin.z0 + 6.0 * 2.54 * 12)
north_shield_rpp = RPP(experiment_rpp.xmin.x0 + 50.0,
                       experiment_rpp.xmax.x0 - 50.0,
                       experiment_rpp.ymax.y0 - 16 * 2.54,
                       experiment_rpp.ymax.y0,
                       experiment_rpp.zmin.z0,
                       experiment_rpp.zmax.z0 + 6.0 * 2.54 * 12)
center_shield_rpp = RPP(west_shield_rpp.xmax.x0 + 5.0 * 2.54 * 12,
                        west_shield_rpp.xmax.x0 + 6.0 * 2.54 * 12,
                        west_shield_rpp.ymin.y0,
                        west_shield_rpp.ymax.y0,
                        west_shield_rpp.zmin.z0,
                        west_shield_rpp.zmax.z0)
east_shield_rpp = RPP(center_shield_rpp.xmax.x0 + 5.0 * 2.54 * 12,
                      center_shield_rpp.xmax.x0 + 6.0 * 2.54 * 12,
                      west_shield_rpp.ymin.y0,
                      west_shield_rpp.ymax.y0,
                      west_shield_rpp.zmin.z0,
                      west_shield_rpp.zmax.z0)
ceiling_plate_rpp = RPP(west_shield_rpp.xmin.x0,
                        east_shield_rpp.xmax.x0,
                        west_shield_rpp.ymin.y0,
                        west_shield_rpp.ymax.y0,
                        west_shield_rpp.zmax.z0,
                        west_shield_rpp.zmax.z0 + 1 * 2.54)
ceiling_shield_rpp = RPP(west_shield_rpp.xmin.x0,
                        east_shield_rpp.xmax.x0,
                        west_shield_rpp.ymin.y0,
                        west_shield_rpp.ymax.y0,
                        ceiling_plate_rpp.zmax.z0,
                        ceiling_plate_rpp.zmax.z0 + 12 * 2.54)

if source_side == 'west':
    libra_center_coord = [(west_shield_rpp.xmax.x0 + center_shield_rpp.xmin.x0)/2,
                          (west_shield_rpp.ymin.y0 + west_shield_rpp.ymax.y0)/2,
                          50]
elif source_side == 'east':
    libra_center_coord = [(center_shield_rpp.xmax.x0 + east_shield_rpp.xmin.x0)/2,
                          (east_shield_rpp.ymin.y0 + east_shield_rpp.ymax.y0)/2,
                          50]

experiment_air_region = -experiment_rpp & +west_shield_rpp & +north_shield_rpp \
     & +center_shield_rpp & +east_shield_rpp & +ceiling_shield_rpp & +ceiling_plate_rpp

west_shield_cell = openmc.Cell(region=-west_shield_rpp, fill=BPE, name='West Shield')
north_shield_cell = openmc.Cell(region=-north_shield_rpp, fill=BPE, name='North Shield')
center_shield_cell = openmc.Cell(region=-center_shield_rpp, fill=BPE, name='Center Shield')
east_shield_cell = openmc.Cell(region=-east_shield_rpp, fill=BPE, name='East Shield')
ceiling_shield_cell = openmc.Cell(region=-ceiling_shield_rpp, fill=BPE, name='Ceiling Shield')
ceiling_plate_cell = openmc.Cell(region=-ceiling_plate_rpp, fill=steel, name='Ceiling Plate')
experiment_air_cell = openmc.Cell(region=experiment_air_region, fill=Air, name='Experiment Air')

universe = openmc.Universe(cells=[
                               west_shield_cell,
                               north_shield_cell,
                               center_shield_cell,
                               east_shield_cell,
                               ceiling_shield_cell,
                               ceiling_plate_cell,
                               experiment_air_cell])

model = build_bates_model(experiment_universe=universe,
                          source_room=source_room,
                          num_particles_per_batch=num_particles_per_batch,
                          libra_center_coord=libra_center_coord)

if __name__ == '__main__':
    curr_dir = os.getcwd()
    os.makedirs(directory, exist_ok=True)
    os.chdir(directory)
    model.export_to_model_xml()
    model.plot_geometry()
    model.run(threads=14)
    os.chdir(curr_dir)


