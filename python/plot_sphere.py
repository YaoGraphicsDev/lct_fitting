import sys
import numpy as np
import pyvista as pv
from charset_normalizer import from_path
from scipy.interpolate import griddata
import argparse


def spherical_to_cartesian(theta, phi):
    """Convert spherical (θ, φ) to Cartesian (x, y, z). θ ∈ [0, π], φ ∈ [0, 2π]"""
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return x, y, z


def load_data(filename):
    """Load [theta_fraction, phi_fraction, value] from a file with auto-detected encoding."""
    result = from_path(filename).best()
    if result is None:
        raise ValueError(f"Could not detect encoding of file: {filename}")

    text = str(result)
    lines = text.splitlines()

    data = []
    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        try:
            parts = [float(x.strip()) for x in line.split(",")]
            if len(parts) == 3:
                data.append(parts)
        except Exception as e:
            print(f"Skipping line: {line} -- {e}")
            continue

    return np.array(data)


def make_colored_sphere(data):
    """Given theta/phi/value data, return a PyVista sphere mesh with interpolated colors."""
    # Convert fractions to radians
    theta = data[:, 0] * np.pi        # [0, π]
    phi = data[:, 1] * np.pi          # [0, 2π]
    values = data[:, 2]

    # Convert data points to Cartesian (unit vectors)
    x_data, y_data, z_data = spherical_to_cartesian(theta, phi)
    points_data = np.stack([x_data, y_data, z_data], axis=1)

    # Create high-res sphere
    sphere = pv.Sphere(radius=1.0, theta_resolution=100, phi_resolution=200)
    points_mesh = sphere.points

    # Keep only hemisphere points: z >= 0
    hemisphere_mask = points_mesh[:, 2] >= 0
    hemisphere_points = points_mesh[hemisphere_mask]

    # Spherical coords for hemisphere points
    xh, yh, zh = hemisphere_points[:, 0], hemisphere_points[:, 1], hemisphere_points[:, 2]
    theta_h = np.arccos(zh)
    phi_h = np.mod(np.arctan2(yh, xh), 2 * np.pi)

    # Interpolate scattered values onto hemisphere points
    interp_points = np.column_stack((theta, phi))
    target_points = np.column_stack((theta_h, phi_h))

    values_on_hemisphere = griddata(interp_points, values, target_points, method='linear')

    # Fill NaNs with nearest neighbor
    nan_mask = np.isnan(values_on_hemisphere)
    if np.any(nan_mask):
        values_on_hemisphere[nan_mask] = griddata(
            interp_points, values, target_points[nan_mask], method='nearest')

    # Assign scalars to sphere mesh
    scalars = np.full(sphere.n_points, np.nan)
    scalars[hemisphere_mask] = values_on_hemisphere
    scalars = np.nan_to_num(scalars, nan=np.nanmin(values_on_hemisphere))

    sphere["values"] = scalars
    return sphere


def add_pole_axis(plotter, center, length=2.5, color="red"):
    start = (center[0], center[1], center[2] - length / 2)  # south pole
    end = (center[0], center[1], center[2] + length / 2)    # north pole
    line = pv.Line(start, end)
    plotter.add_mesh(line, color=color, line_width=3)


def parse_view_direction(view_str):
    """Parse view direction string like '1.0,0.0,0.0' into a camera position."""
    try:
        # Parse the comma-separated values
        view_components = [float(x.strip()) for x in view_str.split(',')]
        if len(view_components) != 3:
            raise ValueError("View direction must have exactly 3 components")
        
        # Normalize the view direction vector
        view_vector = np.array(view_components)
        camera_position = view_vector
        
        # Set up camera position (position, focal point, up vector)
        # For simplicity, we'll use a default up vector of (0, 0, 1)
        # unless the view direction is exactly vertical
        if np.allclose(view_vector, [0, 0, 1]) or np.allclose(view_vector, [0, 0, -1]):
            up_vector = [0, 1, 0]  # Use y-axis as up if looking straight up/down
        else:
            up_vector = [0, 0, 1]  # Default up vector is z-axis
        
        return [camera_position.tolist(), [0, 0, 0], up_vector]
    
    except ValueError as e:
        raise ValueError(f"Invalid view direction format: {e}")


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Plot hemispherical distributions')
    parser.add_argument('datafiles', nargs='+', help='Input data file(s)')
    parser.add_argument('-o', '--output', help='Output image file for headless mode')
    parser.add_argument('-v', '--view', help='Place camera at "x,y,z", looking at 0,0,0 point')
    args = parser.parse_args()
    
    if len(args.datafiles) > 2:
        print("Error: Maximum of 2 data files supported")
        sys.exit(1)
    
    # Parse view direction if provided
    if args.view:
        try:
            camera_position = parse_view_direction(args.view)
        except ValueError as e:
            print(f"Error: {e}")
            sys.exit(1)
    else:
        # Default camera position (current hardcoded value)
        camera_position = [(-10, 0, 5), (0, 0, 0), (0, 0, 1)]
    
    file1 = args.datafiles[0]
    data1 = load_data(file1)
    sphere1 = make_colored_sphere(data1)

    # Set up headless mode if output file is specified
    if args.output:
        plotter = pv.Plotter(off_screen=True)
    else:
        plotter = pv.Plotter()

    if len(args.datafiles) == 2:
        file2 = args.datafiles[1]
        data2 = load_data(file2)
        sphere2 = make_colored_sphere(data2)

        # Translate apart
        sphere1.translate([0, -1.5, 0], inplace=True)
        sphere2.translate([0,  1.5, 0], inplace=True)

        plotter.add_mesh(sphere1, scalars="values", cmap="viridis", show_scalar_bar=True)
        plotter.add_mesh(sphere2, scalars="values", cmap="viridis", show_scalar_bar=True)

        # Labels
        label_points = np.array([
            [0, -1.5, 1.7],
            [0,  1.5, 1.7],
        ])
        labels = [file1, file2]
    else:
        # Only one sphere, centered
        plotter.add_mesh(sphere1, scalars="values", cmap="viridis", show_scalar_bar=True)

        label_points = np.array([[0, 0, 1.7]])
        labels = [file1]

    # Transparent ground plane
    plane = pv.Plane(center=(0, 0, 0), direction=(0, 0, 1), i_size=6, j_size=6)
    plotter.add_mesh(plane, color='lightgray', opacity=0.3)

    # Add pole axis
    if len(args.datafiles) == 2:
        add_pole_axis(plotter, [0, -1.5, 0])
        add_pole_axis(plotter, [0,  1.5, 0])
    else:
        add_pole_axis(plotter, [0, 0, 0])

    plotter.add_axes()
    plotter.camera_position = camera_position

    # Orthogonal camera
    plotter.camera.parallel_projection = True
    bottom, top = -2.0, 2.0
    plotter.camera.SetParallelScale((top - bottom) / 2)

    # Labels
    plotter.add_point_labels(
        label_points, labels,
        font_size=14, text_color='#00FF00',
        point_color='#00FF00', shape=None,
        always_visible=True
    )

    if args.output:
        # Save to file in headless mode
        plotter.screenshot(args.output)
        print(f"Plot saved to {args.output}")
    else:
        # Show interactive window
        plotter.show()


if __name__ == "__main__":
    main()