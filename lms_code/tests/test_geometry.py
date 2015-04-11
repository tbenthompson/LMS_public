from math import sqrt
from lms_code.lib.geometry import *
from lms_code.analysis.collect_data import project_gps

def test_segment_plane_intersection():
    pt1 = (0.0, 0.0)
    pt2 = (3.0, 0.0)
    plane_pt = (2.0, 0.0)
    plane_normal = (1.0, 0.0)
    result = intersection_segment_plane(pt1, pt2, plane_pt, plane_normal)
    np.testing.assert_almost_equal(result[0], [2.0, 0.0])

def test_segment_plane_intersection2():
    pt1 = (0.0, 0.0)
    pt2 = (3.0, 0.0)
    plane_pt = (2.0, 2.0)
    plane_normal = (1.0, 1.0)
    result = intersection_segment_plane(pt1, pt2, plane_pt, plane_normal)
    # Does not intersect
    assert(not result)
    pt2 = (2.0, -2.0)
    result = intersection_segment_plane(pt1, pt2, plane_pt, plane_normal)
    # Check if parallel works
    assert(not result)
    pt2 = (4.0, 0.0)
    result = intersection_segment_plane(pt1, pt2, plane_pt, plane_normal)
    np.testing.assert_almost_equal(result[0], [4.0, 0.0])

    pt1 = (4.0, 0.0)
    pt2 = (6.0, 0.0)
    result = intersection_segment_plane(pt1, pt2, plane_pt, plane_normal)
    np.testing.assert_almost_equal(result[0], [4.0, 0.0])

def test_segment_plane_intersection_in_plane():
    pt1 = (3.0, 1.0)
    pt2 = (1.0, 3.0)
    plane_pt = (2.0, 2.0)
    plane_normal = (1.0, 1.0)
    result = intersection_segment_plane(pt1, pt2, plane_pt, plane_normal)
    np.testing.assert_almost_equal(np.array(result), np.array([pt1,pt2]))

def test_intersection_triangle_plane_nowhere():
    pt1 = (5.0, 0.0)
    pt2 = (6.0, 0.0)
    pt3 = (3.0, 20.0)
    plane_pt = (2.0, 2.0)
    plane_normal = (1.0, 1.0)
    result = intersection_triangle_plane([pt1, pt2, pt3], plane_pt, plane_normal)
    np.testing.assert_almost_equal(result, [])

def test_intersection_triangle_plane_one_point():
    pt1 = (4.0, 0.0)
    pt2 = (6.0, 0.0)
    pt3 = (3.0, 20.0)
    plane_pt = (2.0, 2.0)
    plane_normal = (1.0, 1.0)
    result = intersection_triangle_plane([pt1, pt2, pt3], plane_pt, plane_normal)
    np.testing.assert_almost_equal(result[0, :], [4.0, 0.0])

def test_intersection_triangle_plane_two_points():
    pt1 = (3.0, 0.0)
    pt2 = (6.0, 0.0)
    pt3 = (3.0, 20.0)
    plane_pt = (2.0, 2.0)
    plane_normal = (1.0, 1.0)
    result = intersection_triangle_plane([pt1, pt2, pt3], plane_pt, plane_normal)
    np.testing.assert_almost_equal(result, [[3.0, 1.0], [4.0, 0.0]])

def test_intersection_triangle_plane_in_plane():
    pt1 = (3.0, 1.0)
    pt2 = (1.0, 3.0)
    pt3 = (0.0, 4.0)
    plane_pt = (2.0, 2.0)
    plane_normal = (1.0, 1.0)
    intersection = intersection_triangle_plane([pt1, pt2, pt3],
                                        plane_pt, plane_normal)
    np.testing.assert_almost_equal(intersection, np.array([pt3, pt1, pt2]))

def test_normal():
    normal = plane_normal_from_three_points((1.0, 0.0, 0.0),
                                   (2.0, 0.0, 0.0),
                                   (1.0, 1.0, 0.0))
    np.testing.assert_almost_equal(normal, [0.0, 0.0, 1.0])

def test_project_onto_plane():
    pt1 = (3.0, 0.0)
    plane_pt = (2.0, 2.0)
    plane_normal = (1.0, 1.0)
    projection = project_onto_plane(pt1, plane_pt, plane_normal)
    np.testing.assert_almost_equal(projection, (3.5, 0.5))

def test_project2():
    lon = np.random.rand(50)
    lat = np.random.rand(50)
    proj_lon, proj_lat = project_gps(lon, lat, ((0.5, 0.5, 0.0), (0.5, 0.5, 0.0)))
    for lon, lat in zip(proj_lon, proj_lat):
        np.testing.assert_almost_equal(lon + lat, 1.0)
