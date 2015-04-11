from math import ceil
import numpy as np
from collections import namedtuple
import meshpy.triangle as triangle
from copy import copy
from matplotlib import pyplot as plt
from scipy.sparse import dok_matrix, csgraph
from scipy.spatial import cKDTree
import matplotlib.tri as tri

from codim1.post.interior import interior, interior_on_element
from codim1.core.tools import plot_mesh
from codim1.core import Vertex, Element, PolynomialMapping
from codim1.post import evaluate_solution_on_element

Edge = namedtuple('Edge', 'v_indices, is_boundary, meshpy_idx')

def get_kdtree(mesh, xy_scaling):
    pts = []
    p_per_e = 200
    x_hats = np.linspace(0.0, 1.0, p_per_e)
    for e in mesh:
        for x_hat in x_hats:
            pts.append(e.mapping.get_physical_point(x_hat))
    view_pts = map(lambda x: (x[0], x[1]), pts)
    return cKDTree(np.array(view_pts))

# Need to add and subtract 2 so that the values 0 and 1 are not used for
# boundary markers. These are special and reserved by triangle.
def marker_from_e_idx(e_idx):
    return e_idx + 2

def marker_to_e_idx(e_idx):
    return e_idx - 2

class VolumeMesh(object):
    # viewing_region should be like ((x_min, x_max), (y_min, y_max)
    def __init__(self,
                 bem_mesh,
                 fault_mesh,
                 viewing_region,
                 refine_length = -1,
                 refine_area = 1e9,
                 near_edge_factor = 1.0,
                 extra_edges = []):
        self.mesh = bem_mesh
        self.region = viewing_region
        self.refine_length = refine_length
        self.refine_area = refine_area
        self.near_edge_factor = near_edge_factor

        self.xy_scaling = self.calc_scaling()

        self.elements = copy(self.mesh.elements)
        internal_extra_edges = self.add_view(extra_edges)
        self.process_extra_boundaries(internal_extra_edges)

        # Build the meshpy structures
        self.v_mapping = dict()
        self.marker_to_e = dict()
        self.meshpy_es = []
        self.meshpy_vs = []
        self.meshpy_markers = []
        self.es = []
        self.vs = []
        self.collect()

        # KDTree used for nearest neighbor searches when finding refinements
        self.kdtree = get_kdtree(fault_mesh, self.xy_scaling)

        # Calculate the meshpy triangulation
        self.meshpy()

        # Separate the disjoint subregions
        self.calc_subregions()
        self.identify_regions()

    def min_x(self): return self.region[0][0]
    def min_y(self): return self.region[1][0]
    def max_x(self): return self.region[0][1]
    def max_y(self): return self.region[1][1]
    def calc_scaling(self):
        """ Scaling from physical to screen coordinates. """
        return (self.max_y() - self.min_y()) / float(self.max_x() - self.min_x())



    def collect(self):
        """
        Collect vertices and facets for building an interior mesh from the
        out in.
        """
        for e in self.elements:
            factor = self.refine_factor(e.vertex1, e.vertex2)
            added_verts = self.check_add_vertices(e)
            if added_verts is None:
                continue
            self.add_edges_from_e(e, factor)

    def refine_factor(self, v1, v2):
        # Convert to triangulation coordinates and then calculate edge length
        e_len = np.sqrt((v1.loc[0] * self.xy_scaling -
                         v2.loc[0] * self.xy_scaling) ** 2 +
                        (v1.loc[1] - v2.loc[1]) ** 2)
        # Refine factor is 1 if self.refine_length < e_len or the ratio
        # other (integral!)
        return max(1, ceil(e_len / self.refine_length))

    def add_edges_from_e(self, e, refine_factor):
        # Evenly map the refined vertices along the high order mappings.
        vs_x_hat = np.linspace(0.0, 1.0, refine_factor + 1)
        vs_x = [e.mapping.get_physical_point(x_hat) for x_hat in vs_x_hat]

        # Add points in the case of refine_factor > 1
        vs = [e.vertex1]
        # # Create "Vertex" objects in order to provide an id for each vertex.
        for v_x in vs_x[1:-1]:
            new_v = Vertex(v_x)
            vs.append(new_v)
            other_vert = self.v_mapping.get(new_v.id, None)
            assert(other_vert is None)
            self.add_vertex(new_v)
        vs.append(e.vertex2)

        for i, v in enumerate(vs[:-1]):
            self.add_edge_from_indices([v.id, vs[i + 1].id], e)

    def add_edge_from_indices(self, v_indices, e):
        # Add a volumetric mesh edge from two vertices.
        new_e_idx = len(self.es)
        meshpy_indices = [self.v_mapping[v_id] for v_id in v_indices]
        self.es.append(Edge(meshpy_indices, True, len(self.meshpy_es)))
        self.meshpy_es.append(meshpy_indices)
        self.meshpy_markers.append(marker_from_e_idx(new_e_idx))
        self.marker_to_e[new_e_idx] = e

    def check_add_vertices(self, e):
        # If either of the vertices is in the viewing area, we want the edge
        # If we only take edges that are fully in the viewing area, then
        # intersections with the boundaries will be incomplete.
        either_in = self.in_view(e.vertex1.loc) or self.in_view(e.vertex2.loc)
        if not either_in:
            return None

        vs = [e.vertex1, e.vertex2]
        # Add vertices in case they haven't been added yet (vertices are
        # shared between elements)
        for v in vs:
            if not self.v_mapping.get(v.id, None) is None:
                continue
            self.add_vertex(v)
        return vs

    def add_vertex(self, v):
        #
        new_v_idx = len(self.vs)
        self.vs.append(v)
        self.v_mapping[v.id] = new_v_idx
        self.meshpy_vs.append(v.loc)

    def in_view(self, x):
        return self.min_x() <= x[0] <= self.max_x()\
           and self.min_y() <= x[1] <= self.max_y()

    def create_extra_edge(self, e):
        e.mapping = PolynomialMapping(e)
        self.elements.append(e)

    def process_extra_boundaries(self, boundaries):
        for b in boundaries:
            locs = b[0]
            edges = b[1]
            vs = [Vertex(loc) for loc in locs]
            es = [Element(vs[e[0]], vs[e[1]]) for e in edges]
            for e in es:
                self.create_extra_edge(e)

    def add_view(self, extra_boundaries):
        # Add the rectangular outer boundary of the viewing area.
        view_pts = [
            np.array((self.min_x(), self.min_y())),
            np.array((self.min_x(), self.max_y())),
            np.array((self.max_x(), self.max_y())),
            np.array((self.max_x(), self.min_y()))
        ]
        view_edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
        internal_extra_boundaries = copy(extra_boundaries)
        internal_extra_boundaries.append((view_pts, view_edges))
        return internal_extra_boundaries


    def meshpy(self):
        # Call meshpy and create the delaunay triangulation.

        def centroid(vs):
            return np.mean(vs, axis = 0)

        def distance_to_boundary(pt):
            d, l = self.kdtree.query([pt[0] / self.xy_scaling, pt[1]], k = 1)
            return d

        def refine_func(vertices, area):
            center = centroid(vertices)
            d_bndry = distance_to_boundary(center)
            return bool(area > min(self.refine_area, d_bndry * self.near_edge_factor))

        info = triangle.MeshInfo()

        # Enter triangulation coordinates (so that delaunay angles are
        # reasonable) by multiplying by xy_scaling
        internal_points = map(lambda x: (x[0] * self.xy_scaling, x[1]),
                                   copy(self.meshpy_vs))
        info.set_points(internal_points)
        info.set_facets(self.meshpy_es, facet_markers = self.meshpy_markers)

        mesh = triangle.build(info,
                              refinement_func = refine_func,
                              generate_faces = True)
        self.meshpy = mesh

        self.meshpy_pts = np.array(mesh.points)
        # Exit triangulation coordinates
        self.meshpy_pts[:, 0] /= self.xy_scaling
        self.meshpy_tris = np.array(mesh.elements, dtype = np.int)

    def calc_subregions(self):
        # I calculate the connected components of the viewing region using
        # a graph theoretic approach. This is straightforward since we have
        # edges already. The edges are disconnected at boundaries so that the
        # original boundary mesh disects the viewing region into many areas.
        # This way, we can specify that only the area beneath the surface of
        # the earth should be computed and displayed.
        n_pts = self.meshpy_pts.shape[0]
        connectivity = dok_matrix((n_pts, n_pts))
        for f_idx, f in enumerate(self.meshpy.faces):
            if self.on_boundary(f):
                continue
            connectivity[f[0], f[1]] = 1
        # Connected components are computed using a matrix-based approach in
        # scipy.
        self.n_components, self.components =\
            csgraph.connected_components(connectivity,
                    directed = False,
                    return_labels = True)

    def identify_regions(self):
        self.components = list(self.components)
        self.regions = []
        for r in self.components:
            if self.components.count(r) <= 1:
                continue
            if r in self.regions:
                continue
            self.regions.append(r)
        # TODO: I have a bunch of regions. These are numbered in some unknown
        # fashion. I need to replace all the regions that only have one
        # member with a -1 and then number the remaining regions in ascending
        # order.
        # min_region = min(self.regions)
        # self.regions = [r - min_region for r in self.regions]
        # self.components = [map(lambda c: c - min_region, self.components)

    def on_boundary(self, f):
        # Boundary markers are all greater than 2.
        for i in range(len(f)):
            marker = self.meshpy.point_markers[f[i]]
            if marker >= 2:
                return True
        return False

    def get_evaluator(self):
        ie = InteriorEvaluator(self.meshpy_pts,
                               self.components,
                               self.n_components,
                               self.meshpy_tris,
                               self.regions,
                               self.mesh,
                               self.meshpy.point_markers,
                               self.marker_to_e)
        return ie

class InteriorEvaluator(object):
    def __init__(self, meshpy_pts, components, n_components,
            meshpy_tris, regions, mesh, point_markers, marker_to_e):
        self.meshpy_pts = np.array(meshpy_pts)
        self.components = components
        self.n_components = n_components
        self.meshpy_tris = np.array(meshpy_tris)
        self.regions = regions
        self.mesh = mesh
        self.point_markers = np.array(point_markers)
        self.marker_to_e = marker_to_e

    def viz_vertex_labels(self):
        for i in range(self.meshpy_pts.shape[0]):
            x = self.meshpy_pts[i, 0]
            y = self.meshpy_pts[i, 1]
            label_x_loc = 25
            label_y_loc = 25
            plt.annotate(self.components[i], xy = (x, y),
                         xytext = (label_x_loc, label_y_loc),
                         textcoords = 'offset points',
                         ha = 'right',
                         va = 'bottom',
                         bbox = dict(boxstyle = 'round, pad=0.5',
                                     fc = 'yellow',
                                     alpha = 0.5),
                         arrowprops = dict(arrowstyle = '->',
                                           connectionstyle = 'arc3,rad=0'))

    def viz_mesh(self, selected = False):
        plot_tris = self.meshpy_tris
        if selected:
            plot_tris = self.selected_tris
        plt.triplot(self.meshpy_pts[:, 0], self.meshpy_pts[:, 1], plot_tris)
        # for e in self.es:
        #     pt1 = self.vs[e.v_indices[0]].loc
        #     pt2 = self.vs[e.v_indices[1]].loc
        #     plt.plot([pt1[0], pt2[0]], [pt1[1], pt2[1]], 'k-', linewidth = 6)
        # self.viz_vertex_labels()
        if selected is False:
            for r in self.regions:
                loc = self.region_label_loc(r)
                plt.text(loc[0], loc[1], r, fontsize = 24,
                         bbox=dict(facecolor='red', alpha=0.5))
        plt.show(block = False)

    def in_component(self, tri, comp):
        for i in range(len(tri)):
            p_comp = self.components[tri[i]]
            if p_comp == comp:
                return True
        return False

    def region_label_loc(self, r):
        # Here, I just use the location of the first vertex.
        # I should use some median or mean location. It'd be a bit nicer.
        pos = np.zeros(2)
        n = 0
        for pt_idx, c in enumerate(self.components):
            if c == r:
                pos += self.meshpy_pts[pt_idx]
                n += 1
        pos /= n
        return pos

    def choose_subregion(self, which_component, subset_num, total_subsets):
        if not (0 <= which_component <= self.n_components):
            raise Exception("for choose_subregion, which_component must be" +
                    " a valid component of the interior triangulation")
        #TODO: This ignores the triangles in the corner of a region.
        self.selected_tris = []
        tris_viewed = 0
        for t in self.meshpy_tris:
            if not self.in_component(t, which_component):
                continue
            if tris_viewed % total_subsets == subset_num:
                self.selected_tris.append(t)
            tris_viewed += 1
        self.selected_tris = np.array(self.selected_tris)

    def eval_non_bdry(self, vertex, qs, ek, eval_normal = [0.0, 0.0],
                      type = "displacement"):
        eval = interior(self.mesh, vertex, eval_normal, ek,
                        "soln", type, quad_strategy = qs)
        return eval

    def eval_bdry(self, qs, ek, e, vertex, eval_normal, type, which_side):
        # USES the QBX quadrature algorithm of Klockner 2013
        x0 = e.vertex1.loc[0]
        x1 = e.vertex2.loc[0]
        x_hat = (vertex[0] - x0) / (x1 - x0)
        el_normal = e.mapping.get_normal(x_hat)
        n_dist = 5 * e.length / qs.unit_points

        side_mult = 1
        if type == 'displacement':
            if which_side == "negative":
                side_mult = -1
            if which_side == "positive" and e.bc.type == "traction":
                side_mult = -1
        if type == 'traction':
            side_mult = -1

        c_pt = side_mult * el_normal * n_dist + vertex

        eps = 5e-4
        up_pt = c_pt + side_mult * el_normal * eps
        down_pt = c_pt - side_mult * el_normal * eps

        c_val = self.eval_non_bdry(c_pt, qs, ek, eval_normal, type)
        up_val = self.eval_non_bdry(up_pt, qs, ek, eval_normal, type)
        down_val = self.eval_non_bdry(down_pt, qs, ek, eval_normal, type)

        deriv = (up_val - down_val) / (2 * eps)
        deriv2 = (up_val - 2 * c_val + down_val) / (eps ** 2)
        result = c_val - deriv * (n_dist) - deriv2 * (n_dist ** 2) / 2.0
        return result

    def eval_soln(self, qs, ek, soln_coeffs, which_side, field_normal, field_type):
        soln = dict()
        how_many = 0
        for t in self.selected_tris:
            how_many += 1
            print how_many
            for pt_idx in t:
                vertex = self.meshpy_pts[pt_idx, :]
                if pt_idx in soln.keys():
                    continue
                marker = self.point_markers[pt_idx]
                codim1_e = 0
                if marker >= 2:
                    pre_meshpy_idx = marker_to_e_idx(marker)
                    codim1_e = self.marker_to_e[pre_meshpy_idx]
                if codim1_e != 0 and type(codim1_e.basis) != str:
                    soln[pt_idx] = self.eval_bdry(qs, ek, codim1_e, vertex,
                                             field_normal, field_type, which_side)
                else:
                    soln[pt_idx] = self.eval_non_bdry(vertex, qs, ek,
                                                      field_normal, field_type)
        max_pt_idx = self.meshpy_pts.shape[0]
        return soln
