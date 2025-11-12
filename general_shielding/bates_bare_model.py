import openmc
from openmc.model import RectangularParallelepiped as RPP
import numpy as np

lead = openmc.Material(name="Lead")
lead.add_element('Pb', 1.0,"ao")
lead.set_density("g/cm3", 11.34)

# Name: Borated Polyethylene (5% B in via B4C additive)
# Density: 0.95 g/cm3
# Reference: PNNL Report 15870 (Rev. 1) but revised to make it 5 wt.% B
# Describes: General purpose neutron shielding
BPE = openmc.Material(name='BPE') 
BPE.set_density('g/cm3', 0.95) 
BPE.add_nuclide('H1', 0.1345, 'wo') 
BPE.add_element('B', 0.0500, 'wo') 
BPE.add_element('C', 0.8155, 'wo')

# Portland Concrete from PNNL Materials Compendium (PNNL-15870 Rev2)
PortlandConc = openmc.Material(name="PortlandConc")
PortlandConc.add_element("H", 1e-2, 'wo')
PortlandConc.add_element("C", 0.001000 , 'wo')
PortlandConc.add_element("O", 0.529107, 'wo')
PortlandConc.add_element("Na", 1.6000e-02, 'wo')
PortlandConc.add_element("Mg", 2.0000e-03 , 'wo')
PortlandConc.add_element("Al", 3.3872e-02, 'wo')
PortlandConc.add_element("Si", 3.37021e-01 , 'wo')
PortlandConc.add_element("K", 1.3000e-02, 'wo')
PortlandConc.add_element("Ca", 4.4000e-02 , 'wo')
PortlandConc.add_element("Fe", 1.4000e-02, 'wo')
PortlandConc.set_density('g/cm3', 2.3)

## Steel, Carbon from PNNL Materials Compendium
steel = openmc.Material(name='steel')
steel.set_density('g/cc', 7.82)
steel.add_element('C',  0.0050, 'wo')
steel.add_element('Fe', 0.9950, 'wo')

Air = openmc.Material(name="Air")
Air.add_element("C", 0.00012399 , 'wo')
Air.add_element('N', 0.75527, 'wo')
Air.add_element('O', 0.23178, 'wo')
Air.add_element('Ar', 0.012827, 'wo')
Air.set_density('g/cm3', 0.0012)

# UPDATE WITH ACTUAL DIMENSIONS


def build_bates_model(experiment_universe=None,
                      source_room='warehouse', 
                      run_sim=False,
                      num_particles_per_batch=1e6,
                      libra_center_coord=None,
                      do_translation=True):
    # Coordinates of the center of the bottom of libra tank
    # with the origin at the inside southwest corner of the Bates Lab OC19D

    if source_room == 'warehouse' and libra_center_coord is None:
        libra_center_coord = [80*12*2.54, 30*12*2.54, 50]
    elif source_room == 'cask' and libra_center_coord is None:
        libra_center_coord = [1200, 2250, 50]

    ################### Materials ########################

    # Soil material taken from PNNL Materials Compendium for Earth, U.S. Average
    Soil = openmc.Material(name='Soil')
    Soil.set_density('g/cm3', 1.52)
    Soil.add_element('O', 0.670604, percent_type='ao')
    Soil.add_element('Na', 0.005578, percent_type='ao')
    Soil.add_element('Mg', 0.011432 , percent_type='ao')
    Soil.add_element('Al', 0.053073, percent_type='ao')
    Soil.add_element('Si', 0.201665, percent_type='ao')
    Soil.add_element('K', 0.007653, percent_type='ao')
    Soil.add_element('Ca', 0.026664, percent_type='ao')
    Soil.add_element('Ti', 0.002009, percent_type='ao')
    Soil.add_element('Mn', 0.000272, percent_type='ao')
    Soil.add_element('Fe', 0.021050, percent_type='ao')


    materials = openmc.Materials([steel, PortlandConc, Air, Soil, BPE, lead])

    ################### Geometry #########################

    floor_th = 24*2.54
    outer_wall_th = 0.2
    floor_to_ceil = 30*12*2.54
    ft2cm = 12*2.54
    soil_th = 36*2.54

    corner = {}
    corner['sw'] = [0, 0]
    corner['se'] = [119.7*ft2cm, 0]
    corner['nw'] = [0, 59.6*ft2cm]
    corner['ne'] = [corner['se'][0], corner['nw'][1]]

    ### Surfaces ###

    wall_rpps = {}

    wall_rpps['south'] = RPP(corner['sw'][0] - outer_wall_th, corner['se'][0] + outer_wall_th,
                        corner['sw'][1] - outer_wall_th, corner['se'][1],
                        0, floor_to_ceil)

    wall_rpps['east'] = RPP(corner['se'][0], corner['se'][0] + outer_wall_th,
                        corner['se'][1], corner['ne'][1],
                        0, floor_to_ceil)

    wall_rpps['north'] = RPP(corner['nw'][0] - outer_wall_th, corner['ne'][0] + outer_wall_th,
                        corner['nw'][1], corner['nw'][1] + outer_wall_th,
                        0, floor_to_ceil)

    wall_rpps['west'] = RPP(corner['sw'][0] - outer_wall_th, corner['sw'][0],
                        corner['sw'][1], corner['nw'][1],
                        0, floor_to_ceil)

    floor_rpp = RPP(corner['sw'][0] - outer_wall_th, corner['se'][0] + outer_wall_th,
                    corner['sw'][1] - outer_wall_th, corner['nw'][1] + outer_wall_th,
                    -floor_th, 0)

    ceil_rpp = RPP(floor_rpp.xmin.x0, floor_rpp.xmax.x0,
                floor_rpp.ymin.y0, floor_rpp.ymax.y0,
                floor_to_ceil, floor_to_ceil + outer_wall_th)

    building_air_rpp = RPP(corner['sw'][0], corner['se'][0],
                        corner['sw'][1], corner['nw'][1],
                        0, floor_to_ceil)

    control_outer_rpp = RPP(corner['nw'][0] + 26*ft2cm, corner['nw'][0] + (26 + 16.2)*ft2cm,
                        corner['nw'][1] - (1.1 + 10.3)*ft2cm, corner['nw'][1] - (1.1)*ft2cm,
                        0, 8*ft2cm)
    control_inner_rpp = RPP(control_outer_rpp.xmin.x0 + outer_wall_th, control_outer_rpp.xmax.x0 - outer_wall_th,
                            control_outer_rpp.ymin.y0 + outer_wall_th, control_outer_rpp.ymax.y0 - outer_wall_th,
                            control_outer_rpp.zmin.z0, control_outer_rpp.zmax.z0 - outer_wall_th)
    control_south_concrete_rpp = RPP(control_outer_rpp.xmin.x0, control_outer_rpp.xmin.x0 + 16*ft2cm,
                                     control_outer_rpp.ymin.y0 - 2*ft2cm, control_outer_rpp.ymin.y0,
                                     control_outer_rpp.zmin.z0, control_outer_rpp.zmin.z0 + 8*ft2cm)
    control_east_concrete_rpp = RPP(control_south_concrete_rpp.xmax.x0, control_south_concrete_rpp.xmax.x0 + 2*ft2cm,
                                    control_south_concrete_rpp.ymin.y0 + 9.5*2.54, 
                                    control_south_concrete_rpp.ymin.y0 + 9.5*2.54 + 8*ft2cm,
                                    control_south_concrete_rpp.zmin.z0, control_south_concrete_rpp.zmax.z0)
    

    cask_ent_west_rpp = RPP(corner['nw'][0] + 20.6*ft2cm, corner['nw'][0] + 20.6*ft2cm + 5*2.54,
                            wall_rpps['north'].ymax.y0, wall_rpps['north'].ymax.y0 + 6*ft2cm,
                            0, 10*ft2cm)

    cask_ent_east_rpp = RPP(cask_ent_west_rpp.xmax.x0 + 4*ft2cm, cask_ent_west_rpp.xmax.x0 + 4*ft2cm + 5*2.54,
                            wall_rpps['north'].ymax.y0, wall_rpps['north'].ymax.y0 + 6*ft2cm,
                            0, 10*ft2cm)
    cask_ent_ceil_rpp = RPP(cask_ent_west_rpp.xmin.x0, cask_ent_east_rpp.xmax.x0,
                            cask_ent_west_rpp.ymin.y0, cask_ent_west_rpp.ymax.y0,
                            cask_ent_west_rpp.zmax.z0, cask_ent_west_rpp.zmax.z0 + 5*2.54)
    cask_ent_floor_rpp = RPP(cask_ent_west_rpp.xmin.x0, cask_ent_east_rpp.xmax.x0,
                            cask_ent_west_rpp.ymin.y0, cask_ent_west_rpp.ymax.y0,
                            cask_ent_west_rpp.zmin.z0 - floor_th, cask_ent_west_rpp.zmin.z0)


    cask_west_wall_rpp = RPP(cask_ent_west_rpp.xmax.x0 + 5*2.54 - 77, cask_ent_west_rpp.xmax.x0 + 5*2.54,
                            cask_ent_west_rpp.ymax.y0, cask_ent_west_rpp.ymax.y0 + 411,
                            0, 240)
    cask_north_wall_rpp = RPP(cask_west_wall_rpp.xmin.x0, cask_west_wall_rpp.xmin.x0 + 978,
                            cask_west_wall_rpp.ymax.y0, cask_west_wall_rpp.ymax.y0 + 77,
                            cask_west_wall_rpp.zmin.z0, cask_west_wall_rpp.zmax.z0)
    cask_east_wall_1_rpp = RPP(cask_north_wall_rpp.xmax.x0 - 77, cask_north_wall_rpp.xmax.x0,
                            cask_north_wall_rpp.ymin.y0 - 335, cask_north_wall_rpp.ymin.y0,
                            cask_west_wall_rpp.zmin.z0, cask_west_wall_rpp.zmax.z0)
    cask_south_wall_rpp = RPP(cask_ent_east_rpp.xmin.x0 - 5*2.54, cask_east_wall_1_rpp.xmax.x0,
                            cask_ent_east_rpp.ymax.y0, cask_east_wall_1_rpp.ymin.y0,
                            cask_west_wall_rpp.zmin.z0, cask_west_wall_rpp.zmax.z0)
    cask_east_wall_2_rpp = RPP(cask_south_wall_rpp.xmin.x0, cask_south_wall_rpp.xmin.x0 + 76,
                            cask_south_wall_rpp.ymax.y0, cask_south_wall_rpp.ymax.y0 + 213,
                            cask_west_wall_rpp.zmin.z0, cask_west_wall_rpp.zmax.z0)
    cask_ceil_rpp = RPP(cask_west_wall_rpp.xmin.x0, cask_east_wall_1_rpp.xmax.x0,
                        cask_south_wall_rpp.ymin.y0, cask_north_wall_rpp.ymax.y0,
                        cask_west_wall_rpp.zmax.z0, cask_west_wall_rpp.zmax.z0 + 60)
    cask_floor_rpp = RPP(cask_west_wall_rpp.xmin.x0, cask_east_wall_1_rpp.xmax.x0,
                        cask_south_wall_rpp.ymin.y0, cask_north_wall_rpp.ymax.y0,
                        cask_west_wall_rpp.zmin.z0 - floor_th, cask_west_wall_rpp.zmin.z0)
    


    if source_room == 'cask':
        experiment_rpp = RPP(np.ceil(cask_east_wall_2_rpp.xmax.x0), 
                             np.floor(cask_east_wall_1_rpp.xmin.x0),
                             np.ceil(cask_south_wall_rpp.ymax.y0),
                             np.floor(cask_north_wall_rpp.ymin.y0),
                             np.ceil(cask_west_wall_rpp.zmin.z0),
                             np.floor(cask_west_wall_rpp.zmax.z0))
        print("Cask experiment region: ", (-experiment_rpp).bounding_box)
    else:
        experiment_rpp = RPP(libra_center_coord[0] - 300, libra_center_coord[0] + 400,
                            libra_center_coord[1] - 200, libra_center_coord[1] + 200,
                            0.0, 300)

    soil_bot_plane = openmc.ZPlane(-soil_th, boundary_type='vacuum')
    soil_top_plane = openmc.ZPlane(-10)
    boundary_planes = {}
    boundary_planes['south'] = openmc.YPlane(corner['se'][1] - 200, boundary_type='vacuum')
    boundary_planes['east'] = openmc.XPlane(corner['se'][0] + 200, boundary_type='vacuum')
    boundary_planes['north'] = openmc.YPlane(corner['ne'][1] + 1000, boundary_type='vacuum')
    boundary_planes['west'] = openmc.XPlane(corner['nw'][0] - 200, boundary_type='vacuum')
    boundary_planes['top'] = openmc.ZPlane(ceil_rpp.zmax.z0 + 200, boundary_type='vacuum') 

    # libra_cell = openmc.Cell(region=libra_reg, fill=libra_universe, name='Air LIBRA')


    ### Regions ###

    experiment_reg = -experiment_rpp

    building_reg = +floor_rpp.zmin & -ceil_rpp.zmax \
                & +floor_rpp.ymin & -floor_rpp.ymax \
                & +floor_rpp.xmin & -floor_rpp.xmax

    control_room_wall_reg = -control_outer_rpp & +control_inner_rpp
    control_room_south_concrete_reg = -control_south_concrete_rpp
    control_room_east_concrete_reg = -control_east_concrete_rpp
    control_room_air_reg = -control_inner_rpp

    building_air_reg = -building_air_rpp & +control_outer_rpp & ~experiment_reg \
                    & +control_south_concrete_rpp & +control_east_concrete_rpp

    cask_ent_west_reg = -cask_ent_west_rpp
    cask_ent_east_reg = -cask_ent_east_rpp
    cask_ent_ceil_reg = -cask_ent_ceil_rpp
    cask_west_wall_reg = -cask_west_wall_rpp
    cask_north_wall_reg = -cask_north_wall_rpp
    cask_east_wall_1_reg = -cask_east_wall_1_rpp
    cask_south_wall_reg = -cask_south_wall_rpp | -cask_east_wall_2_rpp
    cask_room_ceil_reg = -cask_ceil_rpp
    cask_floor_reg = -cask_ent_floor_rpp | -cask_floor_rpp

    cask_ent_reg = +cask_ent_west_rpp.xmin & -cask_ent_east_rpp.xmax \
                & +cask_ent_west_rpp.ymin & -cask_ent_west_rpp.ymax \
                & +cask_ent_floor_rpp.zmin & -cask_ent_ceil_rpp.zmax

    cask_room_reg = +cask_west_wall_rpp.xmin & -cask_east_wall_1_rpp.xmax \
                & +cask_south_wall_rpp.ymin & -cask_north_wall_rpp.ymax \
                & +cask_floor_rpp.zmin & -cask_ceil_rpp.zmax

    cask_reg = cask_ent_reg | cask_room_reg

    cask_ent_air_reg = +cask_ent_west_rpp.xmax & -cask_ent_east_rpp.xmin \
                    & +cask_ent_west_rpp.ymin & -cask_ent_west_rpp.ymax \
                    & +cask_floor_rpp.zmax & -cask_ent_ceil_rpp.zmin

    cask_room_air_reg_1 = +cask_west_wall_rpp.xmax & -cask_east_wall_2_rpp.xmin \
                        & +cask_west_wall_rpp.ymin & -cask_north_wall_rpp.ymin \
                        & +cask_floor_rpp.zmax & -cask_ceil_rpp.zmin
    cask_room_air_reg_2 = +cask_east_wall_2_rpp.xmin & -cask_east_wall_2_rpp.xmax \
                        & +cask_east_wall_2_rpp.ymax & -cask_north_wall_rpp.ymin \
                        & +cask_floor_rpp.zmax & -cask_ceil_rpp.zmin
    cask_room_air_reg_3 = +cask_east_wall_2_rpp.xmax & -cask_east_wall_1_rpp.xmin \
                        & +cask_south_wall_rpp.ymax & -cask_north_wall_rpp.ymin \
                        & +cask_floor_rpp.zmax & -cask_ceil_rpp.zmin
    cask_room_air_reg = (cask_room_air_reg_1 | cask_room_air_reg_2 | cask_room_air_reg_3) \
                        & ~experiment_reg

    cask_air_reg = cask_ent_air_reg | cask_room_air_reg


    # print(cask_ent_air_reg.bounding_box)
    # print(cask_room_air_reg_1.bounding_box)
    # print(cask_room_air_reg_2.bounding_box)
    # print(cask_room_air_reg_3.bounding_box)
    # print(cask_room_air_reg.bounding_box)


    soil_reg = -soil_top_plane & +soil_bot_plane \
                & +boundary_planes['west'] & -boundary_planes['east'] \
                & +boundary_planes['south'] & -boundary_planes['north'] \
                & ~cask_floor_reg & +floor_rpp

    outside_air_reg = +soil_top_plane & -boundary_planes['top'] \
                    & +boundary_planes['west'] & -boundary_planes['east'] \
                    & +boundary_planes['south'] & -boundary_planes['north'] \
                    & ~building_reg \
                    & ~cask_reg


    ### Cells ###


    cells = []

    wall_cells = {}
    for direction in wall_rpps.keys():
        wall_cells[direction] = openmc.Cell(region=-wall_rpps[direction], 
                                            fill=steel, 
                                            name='{} wall'.format(direction))
        cells += [wall_cells[direction]]

    floor_cell = openmc.Cell(region=-floor_rpp, fill=PortlandConc, name='floor')
    ceil_cell = openmc.Cell(region=-ceil_rpp, fill=steel, name='ceiling')

    control_room_wall_cell = openmc.Cell(region=control_room_wall_reg, fill=steel, name='control rooom walls')
    control_room_air_cell = openmc.Cell(region=control_room_air_reg, fill=Air, name='control room air')
    control_south_concrete_cell = openmc.Cell(region=control_room_south_concrete_reg, fill=PortlandConc, name='control room south concrete wall')
    control_east_concrete_cell = openmc.Cell(region=control_room_east_concrete_reg, fill=PortlandConc, name='control room east concrete wall')


    ## cask cells
    cask_ent_west_cell = openmc.Cell(region=cask_ent_west_reg, fill=PortlandConc, name='cask entrance west wall')
    cask_ent_east_cell = openmc.Cell(region=cask_ent_east_reg, fill=PortlandConc, name='cask entrance east wall')
    cask_ent_ceil_cell = openmc.Cell(region=cask_ent_ceil_reg, fill=PortlandConc, name='cask entrance ceiling')

    cask_floor_cell = openmc.Cell(region=cask_floor_reg, fill=PortlandConc, name='cask floor')
    cask_room_ceil_cell = openmc.Cell(region=cask_room_ceil_reg, fill=PortlandConc, name='cask room ceiling')

    cask_west_wall_cell = openmc.Cell(region=cask_west_wall_reg, fill=PortlandConc, name='cask west wall')
    cask_north_wall_cell = openmc.Cell(region=cask_north_wall_reg, fill=PortlandConc, name='cask north wall')
    cask_east_wall_cell = openmc.Cell(region=cask_east_wall_1_reg, fill=PortlandConc, name='cask east wall')
    cask_south_wall_cell = openmc.Cell(region=cask_south_wall_reg, fill=PortlandConc, name='cask south wall')
    cask_air_cell = openmc.Cell(region=cask_air_reg, fill=Air, name='cask air')

    building_air_cell = openmc.Cell(region=building_air_reg, fill=Air, name='building air')


    soil_cell = openmc.Cell(region=soil_reg, fill=Soil, name='soil')
    outside_air_cell = openmc.Cell(region=outside_air_reg, fill=Air, name='outside air')

    if isinstance(experiment_universe, openmc.Universe):
        experiment_cell = openmc.Cell(region=experiment_reg, 
                                      fill=experiment_universe, 
                                      name='experiment')
        if do_translation:
            experiment_cell.translation = (libra_center_coord[0], libra_center_coord[1], 0)
    else:
        experiment_cell = openmc.Cell(region=experiment_reg, fill=None, name='experiment')

    cells += [floor_cell, ceil_cell, control_room_wall_cell, control_room_air_cell,
            control_south_concrete_cell, control_east_concrete_cell,
            cask_ent_west_cell, cask_ent_east_cell, cask_ent_ceil_cell,
            cask_floor_cell, cask_room_ceil_cell,
            cask_west_wall_cell, cask_north_wall_cell, cask_east_wall_cell,
            cask_south_wall_cell, cask_air_cell,
            building_air_cell, soil_cell, outside_air_cell,
            experiment_cell
            ]

    universe = openmc.Universe(cells=cells)
    geometry = openmc.Geometry(universe)
    geometry.remove_redundant_surfaces()

    print(floor_cell.region.bounding_box)
    ################### Tallies #########################

    mesh = openmc.RegularMesh(mesh_id=1)
    # mesh.dimension = (27, 27, 22)
    mesh.dimension = (78, 62, 20)
    mesh.lower_left = (-100, -100, 0)
    mesh.upper_right = (3800, 3000, 1000)
    mesh_filter = openmc.MeshFilter(mesh)


    dose_n_energies, dose_n_coeffs = openmc.data.dose_coefficients('neutron')
    # Get rid of coefficients for energies above 15 MeV
    n_mask = dose_n_energies<1.6e7
    dose_n_energies = dose_n_energies[n_mask]
    dose_n_coeffs = dose_n_coeffs[n_mask]

    dose_p_energies, dose_p_coeffs = openmc.data.dose_coefficients('photon')
    # Get rid of coefficients for energies above 20 MeV
    p_mask = dose_p_energies <= 2.0e7
    dose_p_energies = dose_p_energies[p_mask]
    dose_p_coeffs = dose_p_coeffs[p_mask]

    dose_factor =   1e-7 * 3600
    n_dose_filter = openmc.EnergyFunctionFilter(dose_n_energies, dose_n_coeffs * dose_factor, interpolation='histogram')
    p_dose_filter = openmc.EnergyFunctionFilter(dose_p_energies, dose_p_coeffs * dose_factor, interpolation='histogram')

    neutron_filter = openmc.ParticleFilter('neutron')
    photon_filter = openmc.ParticleFilter('photon')

    neutron_tally = openmc.Tally(tally_id=1)
    # neutron_tally.filters.append(mesh_filter)
    neutron_tally.filters.append(mesh_filter)
    neutron_tally.filters.append(n_dose_filter)
    neutron_tally.filters.append(neutron_filter)
    neutron_tally.scores = ['flux']

    photon_tally = openmc.Tally(tally_id=2)
    # photon_tally.filters.append(mesh_filter)
    photon_tally.filters.append(mesh_filter)
    photon_tally.filters.append(p_dose_filter)
    photon_tally.filters.append(photon_filter)
    photon_tally.scores = ['flux']

    # t_tally = openmc.Tally(name='tritium tally')
    # salt_filter = openmc.CellFilter([salt_cell])
    # t_tally.filters.append(salt_filter)
    # t_tally.scores = ['(n,Xt)']

    tallies = openmc.Tallies([neutron_tally, photon_tally])
    # tallies = openmc.Tallies([neutron_tally])
    # tallies = openmc.Tallies()

    ################### Settings #########################

    print("Source Point: ", (libra_center_coord[0] + 0.1, libra_center_coord[1] + 0.1, libra_center_coord[2]))

    point = openmc.stats.Point((libra_center_coord[0] + 0.1, 
                                libra_center_coord[1] + 0.1, 
                                libra_center_coord[2]))
    src = openmc.IndependentSource(space=point)
    src.energy = openmc.stats.Discrete([14.1E6], [1.0])
    src.strength = 1.0

    print(src.space.xyz)
    # vol = openmc.VolumeCalculation(domains=libra_quarter_1_cells, samples=int(1e7))
    settings = openmc.Settings()
    settings.run_mode = 'fixed source'
    settings.source = src
    settings.batches = 100
    settings.inactive = 0
    settings.particles = int(num_particles_per_batch)
    settings.verbosity = 7
    settings.statepoint = {'batches':[ 100]}
    # settings.volume_calculations = [vol]
    settings.photon_transport = True
    # settings.photon_transport = False

    plot_xy = openmc.Plot()
    plot_xy.origin = (libra_center_coord[0], libra_center_coord[1], libra_center_coord[2])
    if source_room == 'warehouse':
        plot_xy.width = (5000, 4000)
        plot_xy.pixels = (5000, 4000)
    elif source_room == 'cask':
        plot_xy.width = (3000, 2000)
        plot_xy.pixels = (3000, 2000)
    plot_xy.basis = 'xy'
    plot_xy.color_by = 'material'

    plot_xz = openmc.Plot()
    plot_xz.origin = (libra_center_coord[0], libra_center_coord[1], libra_center_coord[2] + 500)
    if source_room == 'warehouse':
        plot_xz.width = (5000, 3000)
        plot_xz.pixels = (5000, 3000)
    elif source_room == 'cask':
        plot_xz.width = (3000, 3000)
        plot_xz.pixels = (3000, 3000)
    plot_xz.basis = 'xz'
    plot_xz.color_by = 'material'

    for plot in [plot_xy, plot_xz]:
    #     plot.colors = {
    #     floor_cell: 'darkgray',            # floor (concrete)
    #     ceil_cell: 'slategray',            # ceiling (steel)
    #     control_room_wall_cell: 'dimgray', # control room walls (steel)
    #     control_room_air_cell: 'lightblue',
    #     cask_ent_west_cell: 'peru',
    #     cask_ent_east_cell: 'sandybrown',
    #     cask_ent_ceil_cell: 'burlywood',
    #     cask_floor_cell: 'tan',
    #     cask_room_ceil_cell: 'rosybrown',
    #     cask_west_wall_cell: 'sienna',
    #     cask_north_wall_cell: 'darkred',
    #     cask_east_wall_cell: 'firebrick',
    #     cask_south_wall_cell: 'orangered',
    #     cask_air_cell: 'lightcyan',
    #     building_air_cell: 'lavender',
    #     soil_cell: 'saddlebrown',
    #     outside_air_cell: 'white',
    #     experiment_cell: 'red',
    # }
        plot.colors = {
            PortlandConc: 'silver',
            steel: 'dimgray',
            Air: 'lightcyan',
            Soil: 'saddlebrown',
            BPE: 'limegreen',
            lead: 'gold',
        }

    plots = openmc.Plots([plot_xy, plot_xz])

    # Return the model
    model = openmc.Model(materials=materials, geometry=geometry, settings=settings,
                         tallies=tallies, plots=plots)
    return model

# Example usage:
if __name__ == '__main__':
    model = build_bates_model(source_room='cask')
    model.export_to_model_xml()
    model.plot_geometry()
    # model.run() # Uncomment to run simulation
