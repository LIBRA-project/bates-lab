import openmc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os
from pathlib import Path


# source_room = 'cask' # Options: 'warehouse', 'cask'
# source_side = 'east'
# version_string = 'sides_12in_BPE_north_16in_BPE_ceil_12in_BPE_v1'

source_room = 'cask' # Options: 'warehouse', 'cask'
source_side = 'west'
# version_string = 'sides_24in_conc_east_12in_BPE_ceil_12in_BPE_v1'
# version_string = 'sides_24in_conc_no_east_BPE_ceil_18in_BPE_v1'
version_string = 'center_shield_24in'
directory = Path(source_room) / version_string / source_side


source_rate = 7.0e9 # n/s


ft2cm = 12*2.54
corner = {}
corner['sw'] = [0, 0]
corner['se'] = [119.7*ft2cm, 0]
corner['nw'] = [0, 59.6*ft2cm]
corner['ne'] = [corner['se'][0], corner['nw'][1]]
fence_west_x = corner['nw'][0]
fence_east_x = corner['nw'][0] + 70 * ft2cm
fence_north_y = corner['nw'][1] + 60 * ft2cm

with openmc.StatePoint(directory / 'statepoint.100.h5') as sp:
    n_flux_tally = sp.get_tally(name='neutron flux tally')
    n_flux = n_flux_tally.get_reshaped_data(value='mean', expand_dims=True).squeeze()
    n_flux_err = n_flux_tally.get_reshaped_data(value='std_dev', expand_dims=True).squeeze()

    neutron_tally = sp.get_tally(name='neutron dose tally')
    n_dose = neutron_tally.get_reshaped_data(value='mean', expand_dims=True).squeeze()
    n_dose_err = neutron_tally.get_reshaped_data(value='std_dev', expand_dims=True).squeeze()

    photon_tally = sp.get_tally(name='photon dose tally')
    p_dose = photon_tally.get_reshaped_data(value='mean', expand_dims=True).squeeze()
    p_dose_err = photon_tally.get_reshaped_data(value='std_dev', expand_dims=True).squeeze()


# Plotting
mesh_filter = neutron_tally.find_filter(openmc.MeshFilter)
mesh = mesh_filter.mesh
print("Rectilinear Mesh:", isinstance(mesh, openmc.RectilinearMesh))
print("mesh lower_left:", mesh.lower_left)
print("mesh upper_right:", mesh.upper_right)
xmesh, ymesh = np.meshgrid(np.linspace(mesh.lower_left[0], mesh.upper_right[0], mesh.dimension[0]),
                           np.linspace(mesh.lower_left[1], mesh.upper_right[1], mesh.dimension[1]))
n_dose /= mesh.volumes
p_dose /= mesh.volumes
n_dose_err /= mesh.volumes
p_dose_err /= mesh.volumes

flux_mesh_filter = n_flux_tally.find_filter(openmc.MeshFilter)
flux_mesh = flux_mesh_filter.mesh
n_flux /= mesh.volumes
n_flux_err /= mesh.volumes

total_dose = n_dose + p_dose
total_dose_err = np.sqrt(n_dose_err**2 + p_dose_err**2)
total_dose *= source_rate
total_dose_err *= source_rate

max_ind = np.unravel_index(np.argmax(total_dose), total_dose.shape)
max_position = []
print(max_ind)
if isinstance(mesh, openmc.RegularMesh):
    for i,ind in enumerate(max_ind):
        print(ind)
        max_x = mesh.lower_left[i] + mesh.width[i] * float(ind) + mesh.width[i]/2
        max_position.append(max_x)
elif isinstance(mesh, openmc.RectilinearMesh):
    x_max_position = mesh.x_grid[max_ind[0]] + (mesh.x_grid[max_ind[0]+1] - mesh.x_grid[max_ind[0]])/2
    y_max_position = mesh.y_grid[max_ind[1]] + (mesh.y_grid[max_ind[1]+1] - mesh.y_grid[max_ind[1]])/2
    z_max_position = mesh.z_grid[max_ind[2]] + (mesh.z_grid[max_ind[2]+1] - mesh.z_grid[max_ind[2]])/2
    max_position = [x_max_position, y_max_position, z_max_position]


def plot_total_dose(mesh, basis='xy'):

    if source_room == 'warehouse':
        openmc_overall_plot_width = (
            (mesh.upper_right[0] - mesh.lower_left[0])*1.25,
            (mesh.upper_right[1] - mesh.lower_left[1])*1.25,
            (mesh.upper_right[2] - mesh.lower_left[2])*1.25
        )
        plot_lower_left = mesh.lower_left
        plot_upper_right = mesh.upper_right
    else:
        openmc_overall_plot_width = (3000, 4000, 3000)
        plot_lower_left = (max_position[0] - openmc_overall_plot_width[0]/2,
                           max_position[1] - openmc_overall_plot_width[1]/2,
                           max_position[2] - 200)
        plot_upper_right = (max_position[0] + openmc_overall_plot_width[0]/2,
                            max_position[1] + openmc_overall_plot_width[1]/2,
                            max_position[2] + openmc_overall_plot_width[2] - 200)

    if basis=='xy':
        if isinstance(mesh, openmc.RegularMesh):
            xmesh, ymesh = np.meshgrid(np.linspace(mesh.lower_left[0], mesh.upper_right[0], mesh.dimension[0]),
                                       np.linspace(mesh.lower_left[1], mesh.upper_right[1], mesh.dimension[1]))
        elif isinstance(mesh, openmc.RectilinearMesh):
            xmesh, ymesh = np.meshgrid((mesh.x_grid[:-1] + mesh.x_grid[1:]) / 2, (mesh.y_grid[:-1] + mesh.y_grid[1:]) / 2)

        plot_dose = total_dose[:, :, max_ind[2]].T
        plot_dose_err = total_dose_err[:, :, max_ind[2]].T
        plot_flux = n_flux[:, :, max_ind[2]].T
        openmc_plot_width = openmc_overall_plot_width[0:2]
        xlims = (plot_lower_left[0], plot_upper_right[0])
        ylims = (plot_lower_left[1], plot_upper_right[1])
        x_label = 'X (cm)'
        y_label = 'Y (cm)'
    elif basis=='xz':
        if isinstance(mesh, openmc.RegularMesh):
            xmesh, ymesh = np.meshgrid(np.linspace(mesh.lower_left[0], mesh.upper_right[0], mesh.dimension[0]),
                                       np.linspace(mesh.lower_left[2], mesh.upper_right[2], mesh.dimension[2]))
        elif isinstance(mesh, openmc.RectilinearMesh):
            xmesh, ymesh = np.meshgrid((mesh.x_grid[:-1] + mesh.x_grid[1:]) / 2, (mesh.z_grid[:-1] + mesh.z_grid[1:]) / 2)

        plot_dose = total_dose[:, max_ind[1], :].T
        plot_dose_err = total_dose_err[:, max_ind[1], :].T
        plot_flux = n_flux[:, max_ind[1], :].T
        openmc_plot_width = (openmc_overall_plot_width[0], openmc_overall_plot_width[2])
        xlims = (plot_lower_left[0], plot_upper_right[0])
        ylims = (plot_lower_left[2], plot_upper_right[2])
        x_label = 'X (cm)'
        y_label = 'Z (cm)'
    elif basis=='yz':
        if isinstance(mesh, openmc.RegularMesh):
            xmesh, ymesh = np.meshgrid(np.linspace(mesh.lower_left[1], mesh.upper_right[1], mesh.dimension[1]),
                                       np.linspace(mesh.lower_left[2], mesh.upper_right[2], mesh.dimension[2]))
        elif isinstance(mesh, openmc.RectilinearMesh):
            xmesh, ymesh = np.meshgrid((mesh.y_grid[:-1] + mesh.y_grid[1:]) / 2, (mesh.z_grid[:-1] + mesh.z_grid[1:]) / 2)

        plot_dose = total_dose[max_ind[0], :, :].T
        plot_dose_err = total_dose_err[max_ind[0], :, :].T
        plot_flux = n_flux[max_ind[0], :, :].T
        openmc_plot_width = (openmc_overall_plot_width[1], openmc_overall_plot_width[2])
        xlims = (plot_lower_left[1], plot_upper_right[1])
        ylims = (plot_lower_left[2], plot_upper_right[2])
        x_label = 'Y (cm)'
        y_label = 'Z (cm)'
    
    # convert error to relative error
    plot_dose_err = plot_dose_err / plot_dose * 100.0

    im_ratio = plot_dose.shape[0] / plot_dose.shape[1]
    im_ratio = np.diff(ylims) / np.diff(xlims)
    fig, ax = plt.subplots(figsize=(10, 10 * im_ratio[0]))
    fig_err, ax_err = plt.subplots(figsize=(10, 10 * im_ratio[0]))
    c = ax.pcolormesh(xmesh, ymesh, plot_dose, shading='auto', cmap='plasma',
                    norm=LogNorm(vmin=plot_dose[plot_dose>0].min(), vmax=plot_dose.max()))
    c_err = ax_err.pcolormesh(xmesh, ymesh, plot_dose_err, shading='auto', cmap='Reds',
                    norm=LogNorm(vmin=np.nanmin(plot_dose_err[plot_dose_err>0]), 
                                 vmax=np.nanmax(plot_dose_err)))
    ## add contour lines
    clines = ax.contour(xmesh, ymesh, plot_dose, colors='white', linewidths=0.5, 
                        levels=[1e-3, 1e-2, 1e-1, 0.5, 1, 10, 100, 1000, 1e4, 1e5])
    clines_err = ax_err.contour(xmesh, ymesh, plot_dose_err, colors='black', linewidths=0.5, 
                        levels=[1e-3, 1e-2, 1e-1, 1, 5, 10, 100])
    ## add contour labels
    ax.clabel(clines, fmt=lambda x: '{:.0g}'.format(x), colors='white', fontsize=8)
    ax_err.clabel(clines_err, fmt=lambda x: '{:.0g}%'.format(x), colors='black', fontsize=8)
    fig.colorbar(c, ax=ax, label='Dose (mrem/hr)', fraction=0.046*im_ratio, pad=0.04)

    fig_err.colorbar(c_err, ax=ax_err, label='Dose Standard Deviation (mrem/hr)', fraction=0.046*im_ratio, pad=0.04)

    fig_flux, ax_flux = plt.subplots(figsize=(10, 10 * im_ratio[0]))
    c_flux = ax_flux.pcolormesh(xmesh, ymesh, plot_flux, shading='auto', cmap='viridis',
                    norm=LogNorm(vmin=plot_flux[plot_flux>0].min(), vmax=plot_flux.max()))

    ## add openmc geometry outline
    model = openmc.model.Model.from_model_xml(directory / 'model.xml')
    for a in [ax, ax_err, ax_flux]:
        model.geometry.plot(basis=basis,
                            axes=a,
                            origin=max_position,
                            width=openmc_plot_width,
                            pixels=(int(openmc_plot_width[0]), int(openmc_plot_width[1])),
                            outline='only',
                            color_by='cell'
                            )

        # Plot fence line
        if basis=='xy':
            a.plot([fence_west_x, fence_west_x, fence_east_x, fence_east_x],
                    [corner['nw'][1], fence_north_y, fence_north_y, corner['nw'][1]],
                    'k-', linewidth=2)
        if basis=='xz':
            a.plot([fence_west_x, fence_west_x], [0, 2000], 'k-', linewidth=2)
            a.plot([fence_east_x, fence_east_x], [0, 2000], 'k-', linewidth=2)
        if basis=='yz':
            a.plot([fence_north_y, fence_north_y], [0, 2000], 'k-', linewidth=2)

        a.set_xlim(xlims)
        a.set_ylim(ylims)
        a.set_xlabel(x_label)
        a.set_ylabel(y_label)
        a.set_title('Total Dose Distribution (Neutrons + Photons)')
        a.set_aspect('equal')
    
    fig.tight_layout()
    fig_err.tight_layout()
    fig.savefig(directory / f'total_dose_{basis}.png', dpi=300)
    fig_err.savefig(directory / f'total_dose_err_{basis}.png', dpi=300)

    fig_flux.tight_layout()
    fig_flux.savefig(directory / f'neutron_flux_{basis}.png', dpi=300)

    return fig, ax, fig_err, ax_err

fig_xy, ax_xy, fig_err_xy, ax_err_xy = plot_total_dose(mesh, basis='xy')
fig_xz, ax_xz, fig_err_xz, ax_err_xz = plot_total_dose(mesh, basis='xz')
fig_yz, ax_yz, fig_err_yz, ax_err_yz = plot_total_dose(mesh, basis='yz')
plt.show()