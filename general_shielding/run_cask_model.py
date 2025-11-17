import openmc
from openmc.model import RectangularParallelepiped as RPP
import numpy as np
import os
from pathlib import Path
from bates_bare_model import (build_bates_model, 
                              BPE, lead, PortlandConc, steel,
                              Air)
import copy
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

source_room = 'cask'
source_side = 'west'
version_string = 'center_shield_24in'
num_particles_per_batch = 1e5
use_weight_windows = True

directory = Path(source_room) / version_string / source_side

experiment_rpp = RPP(826.0, 1477.0, 2076.0, 2410.0, 0.0, 240.0)


# west_shield_rpp = RPP(experiment_rpp.xmin.x0 + 100, 
#                       experiment_rpp.xmin.x0 + 100 + 12 * 2.54,
#                       experiment_rpp.ymin.y0,
#                       experiment_rpp.ymin.y0 + 4 * 2.54 * 12,
#                       experiment_rpp.zmin.z0,
#                       experiment_rpp.zmin.z0 + 6.0 * 2.54 * 12)
# north_shield_rpp = RPP(experiment_rpp.xmin.x0 + 50.0,
#                        experiment_rpp.xmax.x0 - 50.0,
#                        experiment_rpp.ymax.y0 - 16 * 2.54,
#                        experiment_rpp.ymax.y0,
#                        experiment_rpp.zmin.z0,
#                        experiment_rpp.zmax.z0 + 6.0 * 2.54 * 12)
# center_shield_rpp = RPP(west_shield_rpp.xmax.x0 + 5.0 * 2.54 * 12,
#                         west_shield_rpp.xmax.x0 + 6.0 * 2.54 * 12,
#                         west_shield_rpp.ymin.y0,
#                         west_shield_rpp.ymax.y0,
#                         west_shield_rpp.zmin.z0,
#                         west_shield_rpp.zmax.z0)
# east_shield_rpp = RPP(center_shield_rpp.xmax.x0 + 5.0 * 2.54 * 12,
#                       center_shield_rpp.xmax.x0 + 6.0 * 2.54 * 12,
#                       west_shield_rpp.ymin.y0,
#                       west_shield_rpp.ymax.y0,
#                       west_shield_rpp.zmin.z0,
#                       west_shield_rpp.zmax.z0)
# ceiling_plate_rpp = RPP(west_shield_rpp.xmin.x0,
#                         east_shield_rpp.xmax.x0,
#                         west_shield_rpp.ymin.y0,
#                         west_shield_rpp.ymax.y0,
#                         west_shield_rpp.zmax.z0,
#                         west_shield_rpp.zmax.z0 + 1 * 2.54)
# ceiling_shield_rpp = RPP(west_shield_rpp.xmin.x0,
#                         east_shield_rpp.xmax.x0,
#                         west_shield_rpp.ymin.y0,
#                         west_shield_rpp.ymax.y0,
#                         ceiling_plate_rpp.zmax.z0,
#                         ceiling_plate_rpp.zmax.z0 + 12 * 2.54)
mid_x = (experiment_rpp.xmin.x0 + experiment_rpp.xmax.x0) / 2
center_shield_rpp = RPP(mid_x - 12 * 2.54, mid_x + 12 * 2.54,
                        experiment_rpp.ymin.y0, experiment_rpp.ymax.y0 - 36 * 2.54,
                        experiment_rpp.zmin.z0, experiment_rpp.zmax.z0 - 12.0 * 2.54)

if source_side == 'west':
    # libra_center_coord = [(west_shield_rpp.xmax.x0 + center_shield_rpp.xmin.x0)/2,
    #                       (west_shield_rpp.ymin.y0 + west_shield_rpp.ymax.y0)/2,
    #                       50]
    libra_center_coord = [(experiment_rpp.xmin.x0 + center_shield_rpp.xmin.x0)/2,
                          experiment_rpp.ymin.y0 + 100,
                          100]
elif source_side == 'east':
    # libra_center_coord = [(center_shield_rpp.xmax.x0 + east_shield_rpp.xmin.x0)/2,
    #                       (east_shield_rpp.ymin.y0 + east_shield_rpp.ymax.y0)/2,
    #                       50]
    libra_center_coord = [(center_shield_rpp.xmax.x0 + experiment_rpp.xmax.x0)/2,
                          experiment_rpp.ymin.y0 + 100,
                          100]

# experiment_air_region = -experiment_rpp & +west_shield_rpp & +north_shield_rpp \
#      & +center_shield_rpp & +east_shield_rpp & +ceiling_shield_rpp & +ceiling_plate_rpp

experiment_air_region = -experiment_rpp & +center_shield_rpp

# west_shield_cell = openmc.Cell(region=-west_shield_rpp, fill=BPE, name='West Shield')
# north_shield_cell = openmc.Cell(region=-north_shield_rpp, fill=BPE, name='North Shield')
center_shield_cell = openmc.Cell(region=-center_shield_rpp, fill=BPE, name='Center Shield')
# east_shield_cell = openmc.Cell(region=-east_shield_rpp, fill=BPE, name='East Shield')
# ceiling_shield_cell = openmc.Cell(region=-ceiling_shield_rpp, fill=BPE, name='Ceiling Shield')
# ceiling_plate_cell = openmc.Cell(region=-ceiling_plate_rpp, fill=steel, name='Ceiling Plate')
experiment_air_cell = openmc.Cell(region=experiment_air_region, fill=Air, name='Experiment Air')

universe = openmc.Universe(cells=[
                            #    west_shield_cell,
                            #    north_shield_cell,
                               center_shield_cell,
                            #    east_shield_cell,
                            #    ceiling_shield_cell,
                            #    ceiling_plate_cell,
                               experiment_air_cell])

model, mesh = build_bates_model(experiment_universe=universe,
                          source_room=source_room,
                          num_particles_per_batch=num_particles_per_batch,
                          libra_center_coord=libra_center_coord,
                          do_translation=False,
                          use_weight_windows=True)

if __name__ == '__main__':
    curr_dir = os.getcwd()
    os.makedirs(directory, exist_ok=True)
    os.chdir(directory)

    # Use random ray

    # # From https://fusion-energy.github.io/neutronics-workshop/tasks/task_14_variance_reduction/5_shielded_room_fw_cadis.html
    # rr_model = copy.deepcopy(model)
    # # turn off photons for random ray model
    # rr_model.settings.photon_transport = False

    # # remove dose tallies and add flux tally for weight window generation
    # flux_tally = openmc.Tally(name='flux tally')
    # mesh_filter = openmc.MeshFilter(mesh)
    # flux_tally.filters.append(mesh_filter)
    # particle_filter = openmc.ParticleFilter('neutron')
    # flux_tally.filters.append(particle_filter)
    # flux_tally.scores = ['flux']
    # rr_model.tallies = openmc.Tallies([flux_tally])

    # print("Generating multigroup cross sections for random ray model...")

    # # first generate multigroup cross sections
    # rr_model.convert_to_multigroup(
    #     # I tend to use "stochastic_slab" method here.
    #     # Using the "material_wise" method is more accurate but slower
    #     # In problems where one needs weight windows to solve we don't really want
    #     # the calculation of weight windows to be slow.
    #     # In extreme cases the "material_wise" method could require its own weight windows to solve.
    #     # The "stochastic_slab" method is much faster and works well for most problems.
    #     # more details here https://docs.openmc.org/en/latest/usersguide/random_ray.html#the-easy-way
    #     method="stochastic_slab", 
    #     overwrite_mgxs_library=True,  # overrights the 
    #     nparticles=2000 # this is the default but can be adjusted upward to improve the fidelity of the generated cross section library
    # )

    # # Convert to random ray model
    # rr_model.convert_to_random_ray()
    # rr_model.settings.random_ray['source_region_meshes'] = [(mesh, [rr_model.geometry.root_universe])]
    # rr_model.settings.weight_window_generators = openmc.WeightWindowGenerator(
    #     method='fw_cadis',
    #     mesh=mesh,
    # )
    # print("Running random ray model to generate weight windows...")
    # random_ray_wwg_statepoint = rr_model.run()

    # # plot the weight windows
    # weight_windows = openmc.hdf5_to_wws('weight_windows.h5')
    # print("Max weight window lower bound:", np.max(weight_windows[0].lower_ww_bounds))
    # print("Min weight window lower bound above 0:", np.min(weight_windows[0].lower_ww_bounds[np.nonzero(weight_windows[0].lower_ww_bounds)]))
    # print("mean weight window lower bound:", np.mean(weight_windows[0].lower_ww_bounds))
    # source_z_ind = int(np.floor((libra_center_coord[2] - mesh.lower_left[2]) / mesh.width[2]))
    # ax1 = plt.subplot()
    # im = ax1.imshow(
    #     weight_windows[0].lower_ww_bounds.squeeze()[:, :, source_z_ind].T,
    #     origin='lower',
    #     extent=mesh.bounding_box.extent['xy'],
    #     # norm=LogNorm(vmin=1e-6, vmax=1e1)
    # )

    # plt.colorbar(im, ax=ax1)

    # ax1 = rr_model.plot(
    #     outline='only',
    #     extent=rr_model.bounding_box.extent['xz'],
    #     axes=ax1,  # Use the same axis as ax1\n",
    #     pixels=10_000_000,  #avoids rounded corners on outline
    #     color_by='material',
    # )
    # ax1.set_title("lower_ww_bounds")

    # plt.show()

    # model.settings.weight_windows = weight_windows
    # model.settings.weight_windows_on = True

    if use_weight_windows:
        if os.path.exists('weight_windows.h5'):
            print("Using existing weight windows...")
            weight_windows = openmc.hdf5_to_wws('weight_windows.h5')
            model.settings.weight_windows = weight_windows
            model.settings.weight_windows_on = True
        else:
            print("Generating weight windows using magic method...")
            wwg = openmc.WeightWindowGenerator(
                        method='magic',
                        mesh=mesh,
                        max_realizations=model.settings.batches
                    )
            model.settings.weight_window_generators = wwg

    print("Running main model...")
    model.export_to_model_xml()
    model.plot_geometry()
    model.run(threads=14, geometry_debug=False)
    os.chdir(curr_dir)


