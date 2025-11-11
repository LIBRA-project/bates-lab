import openmc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os
from pathlib import Path


# source_room = 'cask' # Options: 'warehouse', 'cask'
# source_side = 'east'
# version_string = 'sides_12in_BPE_north_16in_BPE_ceil_12in_BPE_v1'

source_room = 'warehouse' # Options: 'warehouse', 'cask'
source_side = 'south'
# version_string = 'sides_24in_conc_east_12in_BPE_ceil_12in_BPE_v1'
version_string = 'sides_24in_conc_no_east_BPE_ceil_18in_BPE_v1'
directory = Path(source_room) / version_string / source_side

if source_room == 'warehouse':
    libra_center_coord = [80*12*2.54, 30*12*2.54, 50]
elif source_room == 'cask':
    libra_center_coord = [1200, 2250, 50]

source_rate = 3.0e9 # n/s

with openmc.StatePoint(directory / 'statepoint.100.h5') as sp:
    neutron_tally = sp.get_tally(id=1)
    n_dose = neutron_tally.get_reshaped_data(value='mean', expand_dims=True).squeeze()
    n_dose_err = neutron_tally.get_reshaped_data(value='std_dev', expand_dims=True).squeeze()

    photon_tally = sp.get_tally(id=2)
    p_dose = photon_tally.get_reshaped_data(value='mean', expand_dims=True).squeeze()
    p_dose_err = photon_tally.get_reshaped_data(value='std_dev', expand_dims=True).squeeze()


# Plotting
mesh_filter = neutron_tally.find_filter(openmc.MeshFilter)
mesh = mesh_filter.mesh
xmesh, ymesh = np.meshgrid(np.linspace(mesh.lower_left[0], mesh.upper_right[0], mesh.dimension[0]),
                           np.linspace(mesh.lower_left[1], mesh.upper_right[1], mesh.dimension[1]))
n_dose /= mesh.volumes
p_dose /= mesh.volumes
n_dose_err /= mesh.volumes
p_dose_err /= mesh.volumes

total_dose = n_dose + p_dose
total_dose_err = np.sqrt(n_dose_err**2 + p_dose_err**2)
total_dose *= source_rate
total_dose_err *= source_rate

max_ind = np.unravel_index(np.argmax(total_dose), total_dose.shape)
max_position = []
print(max_ind)
for i,ind in enumerate(max_ind):
    print(ind)
    max_x = mesh.lower_left[i] + mesh.width[i] * float(ind) + mesh.width[i]/2
    max_position.append(max_x)

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
        openmc_overall_plot_width = (3000, 2000, 2000)
        plot_lower_left = (max_position[0] - openmc_overall_plot_width[0]/2,
                           max_position[1] - openmc_overall_plot_width[1]/2,
                           max_position[2] - 200)
        plot_upper_right = (max_position[0] + openmc_overall_plot_width[0]/2,
                            max_position[1] + openmc_overall_plot_width[1]/2,
                            max_position[2] + openmc_overall_plot_width[2]/2)

    if basis=='xy':
        xmesh, ymesh = np.meshgrid(np.linspace(mesh.lower_left[0], mesh.upper_right[0], mesh.dimension[0]),
                                   np.linspace(mesh.lower_left[1], mesh.upper_right[1], mesh.dimension[1]))
        plot_dose = total_dose[:, :, max_ind[2]].T
        plot_dose_err = total_dose_err[:, :, max_ind[2]].T
        openmc_plot_width = openmc_overall_plot_width[0:2]
        xlims = (plot_lower_left[0], plot_upper_right[0])
        ylims = (plot_lower_left[1], plot_upper_right[1])
        x_label = 'X (cm)'
        y_label = 'Y (cm)'
    elif basis=='xz':
        xmesh, ymesh = np.meshgrid(np.linspace(mesh.lower_left[0], mesh.upper_right[0], mesh.dimension[0]),
                                   np.linspace(mesh.lower_left[2], mesh.upper_right[2], mesh.dimension[2]))
        plot_dose = total_dose[:, max_ind[1], :].T
        plot_dose_err = total_dose_err[:, max_ind[1], :].T
        openmc_plot_width = (openmc_overall_plot_width[0], openmc_overall_plot_width[2])
        xlims = (plot_lower_left[0], plot_upper_right[0])
        ylims = (plot_lower_left[2], plot_upper_right[2])
        x_label = 'X (cm)'
        y_label = 'Z (cm)'
    elif basis=='yz':
        xmesh, ymesh = np.meshgrid(np.linspace(mesh.lower_left[1], mesh.upper_right[1], mesh.dimension[1]),
                                   np.linspace(mesh.lower_left[2], mesh.upper_right[2], mesh.dimension[2]))
        plot_dose = total_dose[max_ind[0], :, :].T
        plot_dose_err = total_dose_err[max_ind[0], :, :].T
        openmc_plot_width = (openmc_overall_plot_width[1], openmc_overall_plot_width[2])
        xlims = (plot_lower_left[1], plot_upper_right[1])
        ylims = (plot_lower_left[2], plot_upper_right[2])
        x_label = 'Y (cm)'
        y_label = 'Z (cm)'

    im_ratio = plot_dose.shape[0] / plot_dose.shape[1]
    im_ratio = np.diff(ylims) / np.diff(xlims)
    fig, ax = plt.subplots(figsize=(10, 10 * im_ratio[0]))
    c = ax.pcolormesh(xmesh, ymesh, plot_dose, shading='auto', cmap='plasma',
                    norm=LogNorm(vmin=plot_dose[plot_dose>0].min(), vmax=plot_dose.max()))
    ## add contour lines
    clines = ax.contour(xmesh, ymesh, plot_dose, colors='white', linewidths=0.5, 
                        levels=[1e-3, 1e-2, 1e-1, 1, 10, 100, 1000, 1e4, 1e5])
    ## add contour labels
    ax.clabel(clines, fmt=lambda x: '{:.0g}'.format(x), colors='white', fontsize=8)
    fig.colorbar(c, ax=ax, label='Dose (mrem/hr)', fraction=0.046*im_ratio, pad=0.04)

    ## add openmc geometry outline
    model = openmc.model.Model.from_model_xml(directory / 'model.xml')
    model.geometry.plot(basis=basis,
                        axes=ax,
                        origin=max_position,
                        width=openmc_plot_width,
                        pixels=(int(openmc_plot_width[0]), int(openmc_plot_width[1])),
                        outline='only',
                        color_by='cell'
                        )
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title('Total Dose Distribution (Neutrons + Photons)')
    ax.set_aspect('equal')

    fig.tight_layout()
    fig.savefig(directory / f'total_dose_{basis}.png', dpi=300)

    return fig, ax

fig_xy, ax_xy = plot_total_dose(mesh, basis='xy')
fig_xz, ax_xz = plot_total_dose(mesh, basis='xz')
fig_yz, ax_yz = plot_total_dose(mesh, basis='yz')
plt.show()