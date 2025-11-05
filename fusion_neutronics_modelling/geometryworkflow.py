import openmc
import openmc.plotter
import paramak
from cad_to_dagmc import CadToDagmc
from pathlib import Path
# import os

# --- Step 0: Setup ---
# Set path for cross-sections.
# This assumes you are running from a location where 'nuclear_data/cross_sections.xml' exists
# as noted in your environment.
openmc.config['cross_sections'] = Path.home() / 'nuclear_data' / 'cross_sections.xml'
print(f"Set openmc.config['cross_sections'] to {openmc.config['cross_sections']}")

# --- NEW: Define and create the output directory ---
output_dir = Path.cwd() / "outputs"
output_dir.mkdir(parents=True, exist_ok=True)
print(f"All outputs will be saved to: {output_dir}")

# We no longer need the os.environ check from the previous version.


# --- Step 1: Create Paramak Model ---
# We will use a 180-degree sector model to show how boundaries work.
print("Creating Paramak model...")

# FIX 1: Define the build lists as variables first
# This allows us to re-use 'radial_build_list' in Step 5
radial_build_list = [
    (paramak.LayerType.GAP, 10),
    (paramak.LayerType.SOLID, 30),
    (paramak.LayerType.SOLID, 50),
    (paramak.LayerType.SOLID, 10),
    (paramak.LayerType.SOLID, 70),
    (paramak.LayerType.SOLID, 20),
    (paramak.LayerType.GAP, 60),
    (paramak.LayerType.PLASMA, 300),
    (paramak.LayerType.GAP, 60),
    (paramak.LayerType.SOLID, 20),
    (paramak.LayerType.SOLID, 110),
    (paramak.LayerType.SOLID, 10),
]

vertical_build_list = [
    (paramak.LayerType.SOLID, 15),
    (paramak.LayerType.SOLID, 80),
    (paramak.LayerType.SOLID, 10),
    (paramak.LayerType.GAP, 50),
    (paramak.LayerType.PLASMA, 700),
    (paramak.LayerType.GAP, 60),
    (paramak.LayerType.SOLID, 10),
    (paramak.LayerType.SOLID, 40),
    (paramak.LayerType.SOLID, 15),
]

my_reactor = paramak.tokamak(
    radial_build=radial_build_list,
    vertical_build=vertical_build_list,
    triangularity=0.55,
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

# Get the material tags from the reactor model.
# Paramak automatically creates these, e.g., 'plasma', 'blanket', etc.
# This is the fix: .names() is the correct method, not .material_tags
reactor_material_tags = my_reactor.names()
print(f"Reactor material tags found: {reactor_material_tags}")

# --- USER REQUEST: Print attributes of my_reactor ---
print(f"\nAttributes of my_reactor object:\n{dir(my_reactor)}\n")


# --- Step 2: Convert CAD to DAGMC (In-Memory) ---
# This uses the method from '2_converting_cad_in_memory.ipynb'
print("Converting CAD to DAGMC...")
converter = CadToDagmc()

# Add the CadQuery assembly from the Paramak object
# FIX 1: Pass 'my_reactor' directly, as it IS the assembly
converter.add_cadquery_object(
    my_reactor,
    material_tags=reactor_material_tags
)

h5m_filename = output_dir / "dagmc_reactor.h5m"

# Export the DAGMC .h5m file
# Using a coarse mesh for speed. For real physics, you'd need finer meshes.
# FIX 2: Removed 'verbose=False' as per user feedback
converter.export_dagmc_h5m_file(
    filename=str(h5m_filename),
    max_mesh_size=20.0,
    min_mesh_size=5.0
)
print(f"Successfully created {h5m_filename}")


# --- Step 3: Define OpenMC Materials (Answering Q3) ---
# We must create one openmc.Material for each tag in 'reactor_material_tags'.
print("Defining OpenMC materials...")
materials_list = []
for tag in reactor_material_tags:
    # Create a new material with the *exact* name as the tag
    mat = openmc.Material(name=tag)

    # In a real model, you'd add specific elements based on the tag
    # The 'paramak.tokamak' builder tags the plasma as 'plasma'
    # and other layers as 'radial_layer_1', 'vertical_layer_1' etc.
    if tag == 'plasma':
        mat.add_element('H', 1.0, 'ao')  # Dummy plasma
        mat.set_density('g/cm3', 1.0e-6)
    else:
        # Default for all other solid layers
        mat.add_element('Fe', 1.0, 'ao')  # Dummy steel
        mat.set_density('g/cm3', 7.8)

    materials_list.append(mat)

# Create the final Materials collection
my_materials = openmc.Materials(materials_list)
# my_materials.export_to_xml()


# --- Step 4: Define OpenMC Geometry ---
# This uses the method from '1_cad_model_simulation_minimal.ipynb'
print("Defining OpenMC geometry...")

# Create the DAGMC Universe
# --- FIX 1 (from error): Set auto_geom_ids=True to resolve cell ID conflict ---
dag_univ = openmc.DAGMCUniverse(
    filename=str(h5m_filename),
    auto_geom_ids=True)  # <-- Convert to string

# FIX 2: Define a large, static bounding box (in cm)
# This is the robust method from '1_cad_model_simulation_minimal.ipynb'
# It avoids the need to dynamically query the (failing) .val() method.
# We pick a size (12m) that will safely contain the ~7.5m radius reactor.
max_boundary = 1200.0

# Because we made a 180-degree sector (cut on the Y=0 plane),
# we must set a reflective boundary on that plane.
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

# Define the bounding cell !!!
bounding_cell = openmc.Cell(
    name='bounding_cell',
    region=+x_min_surf & -x_max_surf & +y_min_surf & -y_max_surf & +z_min_surf & -z_max_surf,
    fill=dag_univ
)

# !!!
my_geometry = openmc.Geometry([bounding_cell])
# my_geometry.export_to_xml()


# --- Step 5: Define OpenMC Settings & Source ---
print("Defining OpenMC settings...")
# Define a point source in the center of the plasma
# We calculate the approx. major radius from the radial_build

# FIX 3: Use the local 'radial_build_list' variable defined in Step 1
radial_build_values = [item[1] for item in radial_build_list]
plasma_index = radial_build_list.index((paramak.LayerType.PLASMA, 300))
plasma_inner_radius = sum(radial_build_values[:plasma_index])
plasma_thickness = radial_build_values[plasma_index]
approx_major_radius = plasma_inner_radius + (plasma_thickness / 2.0)

print(f"Placing source at approximate major radius: {approx_major_radius} cm")
# --- FIX 4 (RuntimeError): Move source slightly off the y=0.0 boundary ---
# A y-coordinate of 0.0 is exactly on the reflective boundary,
# which causes OpenMC to reject all source particles.
# We nudge it 0.01cm in the +y direction to put it unambiguously inside the geometry.
source_location = (approx_major_radius, 0.01, 0.0)

# --- FIX 2 (proactive): Use IndependentSource to avoid FutureWarning ---
my_source = openmc.IndependentSource()
my_source.space = openmc.stats.Point(source_location)
my_source.angle = openmc.stats.Isotropic()
my_source.energy = openmc.stats.Discrete([14.1e6], [1.0])  # 14.1 MeV neutrons

my_settings = openmc.Settings()
my_settings.batches = 10
my_settings.particles = 500
my_settings.run_mode = 'fixed source'
my_settings.source = my_source
# my_settings.export_to_xml()


# --- Step 6: Plotting (Answering Q2) --- !!!
# This answers your question about viewing the geometry in OpenMC
print("Creating geometry plots...")
plot_xy = openmc.Plot()
plot_xy.filename = str(output_dir / 'reactor_plot_xy')
# Use our static boundary for the plot width
plot_xy.width = (max_boundary * 2, max_boundary * 2)
plot_xy.pixels = (1000, 1000)
plot_xy.basis = 'xy'
plot_xy.origin = (0, 0, 0)
plot_xy.color_by = 'material'
plot_xy.colors = {mat: 'gray' for mat in materials_list}  # Simple colors

# Plot from a different basis (side view)
plot_xz = openmc.Plot()
plot_xz.filename = str(output_dir / 'reactor_plot_xz')
# Use our static boundary for the plot width
plot_xz.width = (max_boundary * 2, max_boundary * 2)
plot_xz.pixels = (1000, 1000)
plot_xz.basis = 'xz'
plot_xz.origin = (0, 0, 0)
plot_xz.color_by = 'material'

plots = openmc.Plots([plot_xy, plot_xz])
plots.export_to_xml(path=str(output_dir / 'plots.xml'))

# Run the plotting
openmc.plot_geometry(cwd=str(output_dir))
print(f"Plots created in {output_dir}: reactor_plot_xy.png, reactor_plot_xz.png")


# --- Step 7: Run Simulation ---
# This will only work if you have cross-sections properly set up
# We removed the 'if' check as we are now setting the cross_sections path directly in Step 0
print("Running OpenMC simulation...")
my_model = openmc.Model(
    geometry=my_geometry,
    materials=my_materials,
    settings=my_settings
)
model_path = my_model.run(cwd=str(output_dir))
print(f"Simulation run complete. Results in {model_path}")


print("\nWorkflow complete.")
