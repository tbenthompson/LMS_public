import numpy as np

def intersection_segment_plane(seg_start, seg_end, plane_pt, plane_normal,
                               epsilon = 1e-5):
    """
    Find the intersection of a line segment and a plane. The plane is defined
    by a normal and a point on the plane.
    """
    seg_end = np.array(seg_end)
    seg_start = np.array(seg_start)
    plane_pt = np.array(plane_pt)
    n = np.array(plane_normal)
    n /= np.linalg.norm(n)

    v = seg_end - seg_start
    w = seg_start - plane_pt
    w2 = seg_end - plane_pt

    wdotn = w.dot(n)
    # If the line segment lies within the plane, just return True
    if abs(wdotn) < epsilon and abs(w2.dot(n)) < epsilon:
        outpoints = [seg_start, seg_end]
        return outpoints

    vdotn = v.dot(n)
    # Check if the line is parallel or almost parallel
    if abs(vdotn) < epsilon:
        return []

    # Find the intersection of the line and the plane
    t = -wdotn / vdotn

    # Check if the intersection is within the segment
    if t < 0:
        return []
    elif t > 1:
        return []

    # Return the intersection
    return [t * v + seg_start]

def unique_rows(a):
    b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    _, idx = np.unique(b, return_index=True)

    unique_a = a[idx]
    return unique_a

def intersection_triangle_plane(vertices, plane_pt, plane_normal):
    """
    Find the intersection of a triangle and a plane. The plane is defined
    by a normal and a point on the plane. The intersection between each
    line segment of the triangle and the plane is computed.
    Each row of vertices should represent one vertex of the triangle.

    The output takes the form of zero, one, two, or three points.

    * Zero points implies no intersection.
    * One point implies that only one vertex of the triangle touches the plane
    * Two points imply a normal line segment intersection.
    * Three points imply that the whole triangle lies in the plane.
    """
    vertices = np.array(vertices)
    intersection_pts =\
        intersection_segment_plane(vertices[0, :], vertices[1, :],
                                   plane_pt, plane_normal)
    intersection_pts.extend(
        intersection_segment_plane(vertices[1, :], vertices[2, :],
                                   plane_pt, plane_normal))
    intersection_pts.extend(
        intersection_segment_plane(vertices[2, :], vertices[0, :],
                                   plane_pt, plane_normal))
    intersection_pts = np.array(intersection_pts)
    if intersection_pts.size == 0:
        return intersection_pts
    intersection_pts = unique_rows(intersection_pts)
    return intersection_pts

def plane_normal_from_three_points(pt1, pt2, pt3):
    """
    Computes (pt3 - pt1) X (pt2 - pt1).
    The result is normalized
    """
    pt1 = np.array(pt1)
    pt2 = np.array(pt2)
    pt3 = np.array(pt3)
    vec1 = pt3 - pt1
    vec2 = pt2 - pt1
    normal = np.cross(vec2, vec1)
    normal /= np.linalg.norm(normal)
    return normal

def project_onto_plane(pt, plane_pt, plane_normal):
    """
    Project the pt onto the plane defined by plane_pt, plane_normal
    In other words, I subtract the orthogonal projection onto plane_normal of
    the vector from plane_pt to pt.
    This find the points on the plane that is nearest to pt.
    """
    normalized_normal = np.array(plane_normal) / np.linalg.norm(plane_normal)
    pt = np.array(pt)
    v = pt - np.array(plane_pt)
    normal_project = v.dot(normalized_normal) * normalized_normal
    new_pt = pt - normal_project
    return new_pt, normal_project

def intersection_surface_plane(vertices, triangles, plane_pt, plane_normal):
    """
    vertices = A numpy array, each row of the array is a vertex
    triangles = List of tuples pointing to indices in the vertices array
    plane_pt = A point on the plane
    plane_normal = The normal to the plane

    returns:
    A list of points of intersection
    A list of line segments connecting those points.
    """

    # Intersect with each triangle
    intersection_verts = []
    segments = []
    for t in triangles:
        local_v = [vertices[v, :] for v in t]
        # Find the intersections
        new_pts = intersection_triangle_plane(local_v, plane_pt, plane_normal)
        if new_pts.size == 0:
            continue
        # Add the vertices and segments to our growing lists
        next_vertex = len(intersection_verts)
        segments.append((next_vertex, next_vertex + 1))
        intersection_verts.extend(list(new_pts))
    return intersection_verts, segments
