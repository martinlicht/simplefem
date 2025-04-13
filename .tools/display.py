

# ----------------------- Import necessary ParaView modules ---------------------

import sys
import os
from paraview.simple import *


# ----------------------- Configuration -----------------------------------------

# The first command line argument (after the script name) indicates the file to be loaded.

if len(sys.argv) < 2:
    print("Usage: pvpython <script_name.py> <path_to_your_vtk_file.vtk>")
    sys.exit(1)

for arg in sys.argv:
    print(f"Argument: {arg}")

file_path = sys.argv[1]

other_arguments = sys.argv[2:] if len(sys.argv) > 2 else []

if not os.path.exists(file_path):
    print(f"Error: File not found at '{file_path}'")
    sys.exit(1)

print(f"Loading file: {file_path}")


# --------------------- ParaView Pipeline -------------------------------------------

# Disable automatic camera reset on 'Show'.
# paraview.simple._DisableFirstRenderCameraReset()

# Create a new VTK Legacy Reader for unstructured grids (.vtk), LegacyVTKReader is appropriate.

reader = LegacyVTKReader(FileNames=[file_path])
print("Reader created.")

# Get the active view or create one if none exists.

view = GetActiveViewOrCreate('RenderView')
if not view:
    print("Error: Could not get or create a RenderView.")
    sys.exit(1)

display_props = Show(reader, view)
print("Data shown in view.")

display_props.Representation = 'Surface With Edges'
print("Representating mesh as 'Surface With Edges'.")


# ----------------------- Inspect Data Arrays --------------------------------

# Update the pipeline. This executes the reader partially to read metadata and make data information available.
reader.UpdatePipeline()
print("Pipeline updated to fetch metadata.")


# Get the data information object for the reader's output port 0. This is the correct way to access array information.
# Can specify output port, default is 0
data_info = reader.GetDataInformation(0)


# Get point data array information container.

point_data_info = data_info.GetPointDataInformation()
print("\nAvailable Point Data Arrays:")
num_point_arrays = point_data_info.GetNumberOfArrays()
if num_point_arrays == 0:
    print("  (None)")
else:
    for i in range(num_point_arrays):
        array_info = point_data_info.GetArrayInformation(i)
        if array_info:
            print(f"  - '{array_info.GetName()}' (Components: {array_info.GetNumberOfComponents()})")
        else:
            print(f"  - (Error retrieving point array info at index {i})")


# Get cell data array information container.
cell_data_info = data_info.GetCellDataInformation()
print("\nAvailable Cell Data Arrays:")
num_cell_arrays = cell_data_info.GetNumberOfArrays()
if num_cell_arrays == 0:
    print("  (None)")
else:
    for i in range(num_cell_arrays):
        array_info = cell_data_info.GetArrayInformation(i)
        if array_info:
            print(f"  - '{array_info.GetName()}' (Components: {array_info.GetNumberOfComponents()})")
        else:
            print(f"  - (Error retrieving cell array info at index {i})")


# ----------------------- Set Coloring --------------------------------

# We attempt to color by the first available array.
color_array_name = None
color_attribute_type = None  # Will be 'POINT_DATA' or 'CELL_DATA'

# We prioritize point data scalars/vectors.
# Otherwise, if no point data were chosen, we use cell data scalars/vectors.

if num_point_arrays > 0:
    first_point_array_info = point_data_info.GetArrayInformation(0)
    if first_point_array_info:
        color_array_name = first_point_array_info.GetName()
        color_attribute_type = 'POINT_DATA'

if num_point_arrays == 0 and num_cell_arrays > 0:
    first_cell_array_info = cell_data_info.GetArrayInformation(0)
    if first_cell_array_info:
        color_array_name = first_cell_array_info.GetName()
        color_attribute_type = 'CELL_DATA'


# However, if any data are indicated on commandline, then check whether they are a scalar field. Use that for coloring.
# Since multiple fields can be given, we process them one by one, even if that means overriding options.

for name in other_arguments:

    for index in range(num_point_arrays):
        point_array_info = point_data_info.GetArrayInformation(index)
        if name == point_array_info.GetName() and point_array_info.GetNumberOfComponents() == 1:
            color_array_name = point_array_info.GetName()
            color_attribute_type = 'POINT_DATA'

    for index in range(num_cell_arrays):
        cell_array_info = cell_data_info.GetArrayInformation(index)
        if name == cell_array_info.GetName() and cell_array_info.GetNumberOfComponents() == 1:
            color_array_name = cell_array_info.GetName()
            color_attribute_type = 'CELL_DATA'


# If an array was found, then we apply it ...
if color_array_name and color_attribute_type:

    print(f"\nColoring by: '{color_array_name}' ({color_attribute_type})")

    ColorBy(display_props, (color_attribute_type, color_array_name))

    # Show the scalar bar (legend)
    display_props.SetScalarBarVisibility(view, True)

    # Force an update of the scalar bar to reflect the chosen array
    UpdateScalarBars(view)

else:
    print("\nNo point or cell data arrays found to color by. Displaying solid color.")
    ColorBy(display_props, None)  # Use solid color


# ----------------------------- Display any specific data-----------------------------------

for name in other_arguments:

    for index in range(num_point_arrays):
        point_array_info = point_data_info.GetArrayInformation(index)
        if name == point_array_info.GetName() and point_array_info.GetNumberOfComponents() == 3:
            glyph = Glyph(Input=reader)
            glyph.OrientationArray = ['POINTS', name]
            glyph.ScaleArray = ['POINTS', name]
            glyph.GlyphType = 'Arrow'
            glyph.ScaleFactor = 0.05
            glyph.GlyphMode = 'All Points'

            # Show in view
            view = GetActiveViewOrCreate('RenderView')
            glyph_display = Show(glyph, view)

            # Color by magnitude
            # glyph_display.Representation = 'Surface'
            # glyph_display.ColorArrayName = ['POINTS', name]
            # glyph_display.SetScalarBarVisibility(view, True)

            view.ResetCamera()
            Render()

    for index in range(num_cell_arrays):
        cell_array_info = cell_data_info.GetArrayInformation(index)
        if name == cell_array_info.GetName() and cell_array_info.GetNumberOfComponents() == 3:
            cell_centers = CellCenters(Input=reader)
            glyph = Glyph(Input=cell_centers)
            glyph.OrientationArray = ['POINTS', name]
            glyph.ScaleArray = ['POINTS', name]
            glyph.GlyphType = 'Arrow'
            glyph.ScaleFactor = 0.05
            glyph.GlyphMode = 'All Points'
            # Show in view
            view = GetActiveViewOrCreate('RenderView')
            glyph_display = Show(glyph, view)

            # # Color by magnitude
            # glyph_display.Representation = 'Surface'
            # glyph_display.ColorArrayName = ['POINTS', name]
            # glyph_display.SetScalarBarVisibility(view, True)

            view.ResetCamera()
            Render()


# ----------------------------- Final Rendering -----------------------------------------

view.ResetCamera()
print("Camera reset.")

Render()
print("Render complete.")

print("\nScript finished. If running in pvpython, close the window to exit.")
Interact()
