# %% [markdown]
# 1: Setup and Package Import

# %%
import matplotlib.pyplot as plt
import openmc
import openmc.plotter
import paramak
from cad_to_dagmc import CadToDagmc
from pathlib import Path
from IPython.display import Image, display  # For displaying plots
import neutronics_material_maker as nmm

# --- Step 0: Setup ---
# Set path for cross-sections.
openmc.config['cross_sections'] = Path.home() / 'nuclear_data' / 'cross_sections.xml'
print(f"Set openmc.config['cross_sections'] to {openmc.config['cross_sections']}")

# --- NEW: Define and create the output directory ---
output_dir = Path.cwd().parent / "outputs"
print(output_dir)
output_dir.mkdir(parents=True, exist_ok=True)
print(f"All outputs will be saved to: {output_dir}")

# %% [markdown]
# 2: Reactor Build

# %%


def basic_tokamak(a, t):
    radial_build_list = [
        (paramak.LayerType.GAP, 10),
        (paramak.LayerType.SOLID, 30),
        (paramak.LayerType.SOLID, 50),
        (paramak.LayerType.SOLID, 10),
        (paramak.LayerType.SOLID, (a / 75) * 70),  # breeder
        (paramak.LayerType.SOLID, 20),
        (paramak.LayerType.GAP, 60),
        (paramak.LayerType.PLASMA, 300),
        (paramak.LayerType.GAP, 60),
        (paramak.LayerType.SOLID, 20),
        (paramak.LayerType.SOLID, (a / 75) * 110),  # breeder
        (paramak.LayerType.SOLID, 10),
    ]

    vertical_build_list = [
        (paramak.LayerType.SOLID, 15),
        (paramak.LayerType.SOLID, (a / 75) * 80),  # breeder
        (paramak.LayerType.SOLID, 10),
        (paramak.LayerType.GAP, 50),
        (paramak.LayerType.PLASMA, 700),
        (paramak.LayerType.GAP, 60),
        (paramak.LayerType.SOLID, 10),
        (paramak.LayerType.SOLID, (a / 75) * 40),  # breeder
        (paramak.LayerType.SOLID, 15),
    ]

    my_reactor = paramak.tokamak(
        radial_build=radial_build_list,
        vertical_build=vertical_build_list,
        triangularity=t,
        rotation_angle=180,
        colors={
            'layer_1': (0.4, 0.9, 0.4),  # center column
            'layer_2': (0.6, 0.8, 0.6),  # magnet shield
            'layer_3': (0.1, 0.8, 0.6),  # first wall
            'layer_4': (0.1, 0.1, 0.9),  # breeder
            'layer_5': (0.4, 0.4, 0.8),  # rear wall
            # plasma has 4 numbers as the last number is the transparency
            'plasma': (1., 0.7, 0.8, 0.6),
        },
    )
    reactor_material_tags = my_reactor.names()

    return (my_reactor, reactor_material_tags, radial_build_list)


# %% [markdown]
# 3: CAD to DAGMC Conversion and Mesh Composition

# %%
def converter_code(my_reactor, reactor_material_tags):
    converter = CadToDagmc()

# Add the CadQuery assembly from the Paramak object
    converter.add_cadquery_object(
        my_reactor,
        material_tags=reactor_material_tags
    )

    h5m_filename = output_dir / "dagmc_reactor.h5m"

# Export the DAGMC .h5m file
    converter.export_dagmc_h5m_file(
        filename=str(h5m_filename),
        max_mesh_size=20.0,
        min_mesh_size=5.0
    )

    return (h5m_filename)


# %% [markdown]
# 4: Materials Definition

# %%
def materials_library(reactor_material_tags):
    materials_list = []

    eurofer_mat = nmm.Material.from_library('eurofer')

    for tag in reactor_material_tags:
        # Create a new material with the *exact* name as the tag
        mat = openmc.Material(name=tag)

    # In a real model, you'd add specific elements based on the tag
    # The 'paramak.tokamak' builder tags the plasma as 'plasma'
    # and other layers as 'radial_layer_1', 'vertical_layer_1' etc.
        if tag == 'plasma':
            mat.add_nuclide('H3', 1.0, 'ao')
            mat.add_nuclide('H2', 1.0, 'ao')  # Dummy plasma
            mat.set_density('g/cm3', 1.0e-6)
        elif tag == 'layer_4':
            mat.add_element(
                'Li',
                0.17,
                percent_type='ao',
                enrichment=60,
                enrichment_target='Li6',
                enrichment_type='ao')
            mat.add_element('Pb', 0.83, 'ao')
            mat.set_density('g/cm3', 11)
        # elif tag == 'layer_4':
            # mat.add_element('He', 1, 'ao')
            # mat.set_density('g/cm3', 0.1786)

        else:
            # Eurofer 97, can't get library import to work
            mat.add_element('Fe', 1, 'ao')
            mat.set_density('g/cm3', 7.8)

        materials_list.append(mat)

# Create the final Materials collection
    my_materials = openmc.Materials(materials_list)

    return (my_materials)

# %% [markdown]
# 5: Defining Model Geometry

# %%


def geometry_code(h5m_filename):

    dag_univ = openmc.DAGMCUniverse(
        filename=str(h5m_filename),
        auto_geom_ids=True)

# Define a large, static bounding box (in cm)
    max_boundary = 1200.0

# Set reflective boundary on the Y=0 plane for the 180-degree sector
    min_x = -max_boundary
    max_x = max_boundary
    min_y = 0.0  # <--- This is the cut plane
    max_y = max_boundary
    min_z = -max_boundary
    max_z = max_boundary

# Define the bounding surfaces
    x_min_surf = openmc.XPlane(min_x, boundary_type='vacuum')
    x_max_surf = openmc.XPlane(max_x, boundary_type='vacuum')
    y_min_surf = openmc.YPlane(min_y, boundary_type='reflective')  # <-- REFLECTIVE
    y_max_surf = openmc.YPlane(max_y, boundary_type='vacuum')
    z_min_surf = openmc.ZPlane(min_z, boundary_type='vacuum')
    z_max_surf = openmc.ZPlane(max_z, boundary_type='vacuum')

# Define the bounding cell
    bounding_cell = openmc.Cell(
        name='bounding_cell',
        region=+x_min_surf & -x_max_surf & +y_min_surf & -y_max_surf & +z_min_surf & -z_max_surf,
        fill=dag_univ
    )

    my_geometry = openmc.Geometry([bounding_cell])
    return (my_geometry)

# %% [markdown]
# 6: Source, Settings & Cells

# %%


def settings_code(radial_build_list):

    radial_build_values = [item[1] for item in radial_build_list]
    plasma_index = radial_build_list.index((paramak.LayerType.PLASMA, 300))
    plasma_inner_radius = sum(radial_build_values[:plasma_index])
    plasma_thickness = radial_build_values[plasma_index]
    approx_major_radius = plasma_inner_radius + (plasma_thickness / 2.0)


# Nudge source slightly off the y=0.0 reflective boundary
    source_location = (approx_major_radius, 0.01, 0.0)

    my_source = openmc.IndependentSource()
    my_source.space = openmc.stats.Point(source_location)
    my_source.angle = openmc.stats.Isotropic()
    my_source.energy = openmc.stats.Discrete([14.1e6], [1.0])  # 14.1 MeV neutrons

    my_settings = openmc.Settings()
    my_settings.batches = 10
    my_settings.particles = 500
    my_settings.run_mode = 'fixed source'
    my_settings.source = my_source
    print("Settings defined.")

    tbr_cell_tally = openmc.Tally(name='tbr')
    tbr_cell_tally.scores = ['(n,Xt)']
    tbr_cell_tally.nuclides = ['Li6', 'Li7']

    mesh = openmc.RegularMesh()  # Need to check the qualities of this mesh
    mesh.dimension = [100, 100, 100]
    # x,y,z coordinates start at 0 as this is a sector model. May need to change.
    mesh.lower_left = [0, 0, -350]
    mesh.upper_right = [650, 650, 350]
    mesh_filter = openmc.MeshFilter(mesh)  # creating a mesh

    tbr_mesh_tally = openmc.Tally(name="tbr_on_mesh")
    tbr_mesh_tally.filters = [mesh_filter]
    tbr_mesh_tally.scores = ["(n,Xt)"]

    tallies = openmc.Tallies([tbr_cell_tally, tbr_mesh_tally])
    return (my_settings, tallies)


# %%
avg_thickenss = 75  # blanket thickness in centimeters, standard value of 75cm
triangularity = 0.55  # triangularity, standard value of 0.55


thickness_values = [50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]
tbr_values = []
tbr_errors = []

for i in thickness_values:
    reactor = basic_tokamak(i, triangularity)

    my_model = openmc.Model(
        geometry=geometry_code(converter_code(reactor[0], reactor[1])),
        materials=materials_library(reactor[1]),
        settings=settings_code(reactor[2])[0],
        tallies=settings_code(reactor[2])[1]
    )

    !rm * .h5

    statepoint_file = my_model.run()

    sp = openmc.StatePoint("statepoint.10.h5")

    tbr_cell_tally = sp.get_tally(name='tbr')

    tbr_values.append(tbr_cell_tally.mean.sum())
    tbr_errors.append(tbr_cell_tally.std_dev.sum())


# tbr_mesh = sp.get_tally(name='tbr_on_mesh')


print(f"The reactor has a TBR of {tbr_cell_tally.mean.sum()}")
print(f"Standard deviation on the TBR is {tbr_cell_tally.std_dev.sum()}")


# df = tbr_cell_tally.get_pandas_dataframe()
# df


plt.plot(thickness_values, tbr_values, 'k.')
plt.errorbar(thickness_values, tbr_values, tbr_errors, fmt='k.')

plt.title('TBR with Varying Blanket Thickness')
plt.xlabel('Average Blanket Thickness (cm)')
plt.ylabel('TBR')

plt.show()

# %%
print(tbr_values)
print(tbr_errors)


plt.plot(thickness_values, tbr_values, 'k.')
plt.errorbar(thickness_values, tbr_values, tbr_errors, fmt='k.')

plt.title('TBR with Varying Blanket Thickness')
plt.xlabel('Average Blanket Thickness (cm)')
plt.ylabel('TBR')

plt.show()
