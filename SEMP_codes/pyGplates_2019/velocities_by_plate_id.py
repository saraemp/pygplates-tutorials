
"""
    Copyright (C) 2018 The University of Sydney, Australia
    
    This program is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License, version 2, as published by
    the Free Software Foundation.
    
    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    
    You should have received a copy of the GNU General Public License along
    with this program; if not, write to Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""


####################################################################################
# Calculate velocities by plate ID at points or along tessellated lines over time. #
####################################################################################


from __future__ import print_function
import argparse
import math
import os.path
import pygplates
import sys


# Required pygplates version.
# Use the most recent public release to avoid any prior bugs in pyGPlates.
PYGPLATES_VERSION_REQUIRED = pygplates.Version(18)

# The default threshold sampling distance along subduction zones.
DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES = 0.5
DEFAULT_THRESHOLD_SAMPLING_DISTANCE_KMS = math.radians(DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES) * pygplates.Earth.equatorial_radius_in_kms

DEFAULT_TIME_RANGE_YOUNG_TIME = 0
DEFAULT_TIME_RANGE_OLD_TIME = 200
DEFAULT_TIME_INCREMENT = 1

DEFAULT_VELOCITY_DELTA_TIME = 1


def calculate_velocities(
        # Rotation model or feature collection(s), or list of features, or filename(s)...
        rotation_features_or_model,
        # Geometry feature collection(s), or list of features, or filename(s) or any combination of those...
        reconstructable_features,
        # Threshold sampling distance along subduction zones (in radians)...
        threshold_sampling_distance_radians,
        reconstruction_time,
        velocity_delta_time = 1.0,
        anchor_plate_id = 0):
    """
    Reconstructs geometries at 'reconstruction_time', tessellates all *line* geometries to within
    'threshold_sampling_distance_radians' radians and returns the following velocity-related parameters
    at each tessellates point:
    
    - point longitude
    - point latitude
    - absolute (relative to anchor plate) velocity magnitude (in cm/yr)
    - absolute velocity obliquity angle (angle between line normal vector and absolute velocity vector)
      * note that this is zero if the input geometries are points (not lines)
    - length of arc segment (in degrees) that current point is on
      * note that this is zero if the input geometries are points (not lines)
    - line arc normal azimuth angle (clockwise starting at North, ie, 0 to 360 degrees) at current point
      * note that this is zero if the input geometries are points (not lines)
    - plate ID
    
    The obliquity angles are in the range (-180 180). The range (0, 180) goes clockwise
    (when viewed from above the Earth) from the line normal direction to the velocity vector.
    The range (0, -180) goes counter-clockwise.
    You can change the range (-180, 180) to the range (0, 360) by adding 360 to negative angles.
    
    """
    
    # Turn rotation data into a RotationModel (if not already).
    rotation_model = pygplates.RotationModel(rotation_features_or_model)
    
    # Turn geometry data into a list of features (if not already).
    reconstructable_features = pygplates.FeaturesFunctionArgument(reconstructable_features)
    
    # Reconstruct our reconstructable features to the current 'reconstruction_time'.
    reconstructed_feature_geometries = []
    pygplates.reconstruct(reconstructable_features.get_features(), rotation_model, reconstructed_feature_geometries, reconstruction_time, anchor_plate_id)
    
    # List of tesselated line points and associated velocity parameters for the current 'reconstruction_time'.
    output_data = []
    
    # Iterate over all reconstructed feature geometries.
    for reconstructed_feature_geometry in reconstructed_feature_geometries:
        
        reconstruction_plate_id = reconstructed_feature_geometry.get_feature().get_reconstruction_plate_id()
        
        # Get the rotation of the reconstructed geometry relative to the anchor plate
        # from 'reconstruction_time + velocity_delta_time' to 'reconstruction_time'.
        #
        # The geometries have been reconstructed using the rotation "R(0->t2,A->M)":
        #
        #   reconstructed_geometry = R(0->t2,A->M) * present_day_geometry
        #
        # We can write "R(0->t2,A->M)" in terms of the stage rotation "R(t1->t2,F->M)" as:
        #
        #   R(0->t2,A->M) = R(0->t2,A->F) * R(0->t2,F->M)
        #                 = R(0->t2,A->F) * R(t1->t2,F->M) * R(0->t1,F->M)
        #                 = R(0->t2,A->F) * stage_rotation * R(0->t1,F->M)
        #
        # ...where 't1' is 't+1' and 't2' is 't' (ie, from 't1' to 't2').
        #
        # So to get the *reconstructed* geometry into the stage rotation reference frame
        # we need to rotate it by "inverse[R(0->t2,A->F)]":
        #
        #   reconstructed_geometry = R(0->t2,A->F) * stage_rotation * R(0->t1,F->M) * present_day_geometry
        #   inverse[R(0->t2,A->F)] * reconstructed_geometry = stage_rotation * R(0->t1,F->M) * present_day_geometry
        #
        # Once we've done that we can calculate the velocities of those geometry points
        # using the stage rotation. Then the velocities need to be rotated back from the
        # stage rotation reference frame using the rotation "R(0->t2,A->F)".
        #
        # However, since we're calculating the *absolute* velocity (ie, relative to anchor plate)
        # the above to/from stage rotation frame conversion "R(0->t2,A->F)" is the
        # identity rotation since the fixed plate (F) is the anchor plate (A).
        #
        # This long explanation was obviously unnecessary in our case but it's good to note in case this
        # code ever gets modified to support *convergence* velocities (ie, relative to non-anchor plate).
        equivalent_stage_rotation = rotation_model.get_rotation(
                reconstruction_time,
                reconstruction_plate_id,
                reconstruction_time + velocity_delta_time,
                anchor_plate_id=anchor_plate_id)
        
        calculate_velocities_along_reconstructed_geometry(
            output_data,
            reconstructed_feature_geometry.get_reconstructed_geometry(),
            equivalent_stage_rotation,
            reconstruction_plate_id,
            threshold_sampling_distance_radians,
            velocity_delta_time)
    
    # Return data sorted since it's easier to compare results (when at least lon/lat is sorted).
    return output_data


def calculate_velocities_along_reconstructed_geometry(
        output_data,
        reconstructed_geometry,
        equivalent_stage_rotation,
        reconstruction_plate_id,
        threshold_sampling_distance_radians,
        velocity_delta_time):
    
    # Ensure the reconstructed geometry is tessellated to within the threshold sampling distance.
    # Ignores tessellation if geometry is a point or multipoint.
    try:
        tessellated_reconstructed_geometry = reconstructed_geometry.to_tessellated(threshold_sampling_distance_radians)
        
        # Iterate over the great circle arcs of the tessellated polyline/polygon to get the
        # arc midpoints, lengths and normals.
        # There is an arc between each adjacent pair of points in the polyline/polygon.
        arc_midpoints = []
        arc_lengths = []
        arc_normals = []
        for arc in tessellated_reconstructed_geometry.get_segments():
            if not arc.is_zero_length():
                arc_midpoints.append(arc.get_arc_point(0.5))
                arc_lengths.append(arc.get_arc_length())
                # The normal to the polyline/polygon.
                arc_normals.append(arc.get_great_circle_normal())
        
    except AttributeError:
        arc_midpoints = None
    
    # If either:
    #
    #   - the reconstructed geometry is a point or multipoint, or
    #   - the tessellated polyline/polygon coincided with a point.
    #
    # ...then just use the original points, and set arc lengths and normals to zero.
    if not arc_midpoints:
        
        reconstructed_points = reconstructed_geometry.get_points()
        
        # Calculate the absolute velocities at the reconstructed points.
        absolute_velocity_vectors = pygplates.calculate_velocities(
                reconstructed_points, equivalent_stage_rotation,
                velocity_delta_time, pygplates.VelocityUnits.cms_per_yr)
        
        for point_index, point in enumerate(reconstructed_points):
        
            lat, lon = point.to_lat_lon()
            
            # The data will be output in GMT format (ie, lon first, then lat, etc).
            output_data.append((
                    lon,
                    lat,
                    absolute_velocity_vectors[point_index].get_magnitude(),
                    0.0,  # absolute_obliquity_degrees
                    0.0,  # math.degrees(arc_length)
                    0.0,  # math.degrees(arc_normal_azimuth)
                    reconstruction_plate_id))
        
        return
    
    # The arc normals relative to North (azimuth).
    # Convert global 3D normal vectors to local (magnitude, azimuth, inclination) tuples (one tuple per point).
    arc_local_normals = pygplates.LocalCartesian.convert_from_geocentric_to_magnitude_azimuth_inclination(
            arc_midpoints, arc_normals)
    
    # Calculate the absolute velocities at the arc midpoints.
    absolute_velocity_vectors = pygplates.calculate_velocities(
            arc_midpoints, equivalent_stage_rotation,
            velocity_delta_time, pygplates.VelocityUnits.cms_per_yr)
    
    for arc_index in range(len(arc_midpoints)):
        arc_midpoint = arc_midpoints[arc_index]
        arc_length = arc_lengths[arc_index]
        arc_normal = arc_normals[arc_index]
        arc_normal_azimuth = arc_local_normals[arc_index][1]
        lat, lon = arc_midpoint.to_lat_lon()
        
        # Calculate the absolute rate parameters.
        absolute_velocity_vector = absolute_velocity_vectors[arc_index]
        if absolute_velocity_vector.is_zero_magnitude():
            absolute_velocity_magnitude = 0
            absolute_obliquity_degrees = 0
        else:
            absolute_velocity_magnitude = absolute_velocity_vector.get_magnitude()
            absolute_obliquity_degrees = math.degrees(pygplates.Vector3D.angle_between(
                    absolute_velocity_vector, arc_normal))
            
            # The direction towards which we rotate from the normal in a clockwise fashion.
            clockwise_direction = pygplates.Vector3D.cross(arc_normal, arc_midpoint.to_xyz())
            
            # Anti-clockwise direction has range (0, -180) instead of (0, 180).
            if pygplates.Vector3D.dot(absolute_velocity_vector, clockwise_direction) < 0:
                absolute_obliquity_degrees = -absolute_obliquity_degrees
        
        # The data will be output in GMT format (ie, lon first, then lat, etc).
        output_data.append((
                lon,
                lat,
                absolute_velocity_magnitude,
                absolute_obliquity_degrees,
                math.degrees(arc_length),
                math.degrees(arc_normal_azimuth),
                reconstruction_plate_id))


def write_output_file(output_filename, output_data):
    with open(output_filename, 'w') as output_file:
        for output_line in output_data:
            output_file.write(' '.join(str(item) for item in output_line) + '\n')


def calculate_velocities_over_time(
        output_filename_prefix,
        output_filename_extension,
        rotation_filenames,
        reconstructable_filenames,
        threshold_sampling_distance_radians,
        time_young,
        time_old,
        time_increment,
        velocity_delta_time = 1.0,
        anchor_plate_id = 0):
    
    if time_increment <= 0:
        print('The time increment "{0}" is not positive and non-zero.'.format(time_increment), file=sys.stderr)
        return
    
    if time_young > time_old:
        print('The young time {0} is older (larger) than the old time {1}.'.format(time_young, time_old), file=sys.stderr)
        return
    
    rotation_model = pygplates.RotationModel(rotation_filenames)
    
    # Read/parse the reconstructable features once so we're not doing at each time iteration.
    reconstructable_features = [pygplates.FeatureCollection(reconstructable_filename)
            for reconstructable_filename in reconstructable_filenames]
    
    # Iterate over the time range.
    time = time_young
    while time <= pygplates.GeoTimeInstant(time_old):
        
        print('Time {0}'.format(time))
        
        # Returns a list of tesselated subduction zone points and associated convergence parameters
        # to write to the output file for the current 'time'.
        output_data = calculate_velocities(
                rotation_model,
                reconstructable_features,
                threshold_sampling_distance_radians,
                time,
                velocity_delta_time,
                anchor_plate_id)
        
        if output_data:
            output_filename = '{0}_{1:0.2f}.{2}'.format(output_filename_prefix, time, output_filename_extension)
            write_output_file(output_filename, output_data)

        # Increment the time further into the past.
        time += time_increment
    
    return 0 # Success


if __name__ == '__main__':
    
    # Check the imported pygplates version.
    if not hasattr(pygplates, 'Version') or pygplates.Version.get_imported_version() < PYGPLATES_VERSION_REQUIRED:
        print('{0}: Error - imported pygplates version {1} but version {2} or greater is required'.format(
                os.path.basename(__file__), pygplates.Version.get_imported_version(), PYGPLATES_VERSION_REQUIRED),
            file=sys.stderr)
        sys.exit(1)
    
    
    __description__ = \
    """Calculate velocities by plate ID at points or along tessellated lines over time.
    
    For each time (over a range of times) an output xy file is generated containing reconstructed points
    (with point locations as the first two columns x and y) and the following parameters in subsequent columns:
    
      - absolute (relative to anchor plate) velocity magnitude (in cm/yr)
      - absolute velocity obliquity angle (angle between line normal vector and absolute velocity vector)
        * note that this is zero if the input geometries are points (not lines)
      - length of arc segment (in degrees) that current point is on
        * note that this is zero if the input geometries are points (not lines)
      - line arc normal azimuth angle (clockwise starting at North, ie, 0 to 360 degrees) at current point
        * note that this is zero if the input geometries are points (not lines)
      - plate ID
    
    The obliquity angles are in the range (-180 180). The range (0, 180) goes clockwise (when viewed from above the Earth) from the
    line normal direction to the velocity vector. The range (0, -180) goes counter-clockwise.
    You can change the range (-180, 180) to the range (0, 360) by adding 360 to negative angles.
    
    Each point in the output is the midpoint of a great circle arc between two adjacent tessellated points in a reconstructed polyline/polygon.
    The normal vector used in the obliquity calculations is perpendicular to the great circle arc of each point (arc midpoint).
    
    Each polyline/polygon is sampled at approximately uniform intervals along its length (specified via a threshold sampling distance).
    The sampling along the entire length of a polyline/polygon is not exactly uniform. Each segment is sampled such that the samples have
    a uniform spacing that is less than or equal to the threshold sampling distance. However each segment might have a slightly different
    spacing distance (since segment lengths are not integer multiples of the threshold sampling distance).

    NOTE: Separate the positional and optional arguments with '--' (workaround for bug in argparse module).
    For example...

    python %(prog)s -r rotations.rot -m reconstructable_features.gpml -t 0 200 -i 1 -v 1 -d 0.5 -e xy -- velocity
     """

    # The command-line parser.
    parser = argparse.ArgumentParser(description = __description__, formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-r', '--rotation_filenames', type=str, nargs='+', required=True,
            metavar='rotation_filename', help='One or more rotation files.')
    parser.add_argument('-m', '--reconstructable_filenames', type=str, nargs='+', required=True,
            metavar='reconstructable_filename', help='One or more reconstructable files to generate reconstructed geometries.')
    parser.add_argument('-a', '--anchor', type=int, default=0,
            dest='anchor_plate_id',
            help='Anchor plate id used for reconstructing. Defaults to zero.')
    
    # Can specify only one of '-i', '-l' or '-t'.
    threshold_sampling_distance_group = parser.add_mutually_exclusive_group()
    threshold_sampling_distance_group.add_argument('-d', '--threshold_sampling_distance_degrees', type=float,
            help='Threshold sampling distance along subduction zones (in degrees). '
                'Defaults to {0} degrees.'.format(DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES))
    threshold_sampling_distance_group.add_argument('-k', '--threshold_sampling_distance_kms', type=float,
            help='Threshold sampling distance along polylines/polygons (in Kms). '
                'Defaults to {0:.2f} Kms (which is equivalent to {1} degrees).'.format(
                        DEFAULT_THRESHOLD_SAMPLING_DISTANCE_KMS,
                        DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES))

    parser.add_argument('-t', '--time_range', type=float, nargs=2,
            metavar=('young_time', 'old_time'),
            default=[DEFAULT_TIME_RANGE_YOUNG_TIME, DEFAULT_TIME_RANGE_OLD_TIME],
            help='The time range (in Ma) from young time to old time. '
                'Defaults to {0} -> {1} Ma.'.format(
                    DEFAULT_TIME_RANGE_YOUNG_TIME, DEFAULT_TIME_RANGE_OLD_TIME))
    
    def parse_positive_number(value_string):
        try:
            value = float(value_string)
        except ValueError:
            raise argparse.ArgumentTypeError("%s is not a number" % value_string)
        
        if value <= 0:
            raise argparse.ArgumentTypeError("%g is not a positive number" % value)
        
        return value
    
    parser.add_argument('-i', '--time_increment', type=parse_positive_number,
            default=DEFAULT_TIME_INCREMENT,
            help='The time increment in My. Defaults to {0} My.'.format(DEFAULT_TIME_INCREMENT))
    
    parser.add_argument('-v', '--velocity_delta_time', type=parse_positive_number,
            default=DEFAULT_VELOCITY_DELTA_TIME,
            help='The delta time interval used to calculate velocities in My. '
                'Defaults to {0} My.'.format(DEFAULT_VELOCITY_DELTA_TIME))
    
    parser.add_argument('output_filename_prefix', type=str,
            help='The output filename prefix. An output file is created for each geological time in the sequence where '
                'the filename suffix contains the time and the filename extension.')
    parser.add_argument('-e', '--output_filename_extension', type=str, default='xy',
            help='The output xy filename extension. Defaults to "xy".')
    
    # Parse command-line options.
    args = parser.parse_args()
    
    if args.time_range[0] > args.time_range[1]:
        raise argparse.ArgumentTypeError("First (young) value in time range is greater than second (old) value")
    
    # Determine threshold sampling distance.
    if args.threshold_sampling_distance_degrees:
        threshold_sampling_distance_radians = math.radians(args.threshold_sampling_distance_degrees)
    elif args.threshold_sampling_distance_kms:
        threshold_sampling_distance_radians = args.threshold_sampling_distance_kms / pygplates.Earth.equatorial_radius_in_kms
    else: # default...
        threshold_sampling_distance_radians = math.radians(DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES)
    
    return_code = calculate_velocities_over_time(
            args.output_filename_prefix,
            args.output_filename_extension,
            args.rotation_filenames,
            args.reconstructable_filenames,
            threshold_sampling_distance_radians,
            args.time_range[0],
            args.time_range[1],
            args.time_increment,
            args.velocity_delta_time,
            args.anchor_plate_id)
    if return_code is None:
        sys.exit(1)
        
    sys.exit(0)
