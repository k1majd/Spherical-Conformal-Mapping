import halfedge_mesh
from halfedge_mesh import Vertex, Facet, Halfedge
import math


class ConformalMap(halfedge_mesh.HalfedgeMesh):
    """Conformal mapping super class of halfedge mesh"""

    def __init__(self, mesh_path):
        super().__init__(mesh_path)

    def star_map(self):
        """map the vertices in the star of a vertex"""
        avg_point = [0.0, 0.0, 0.0]
        for vertex in self.vertices:
            avg_point[0] += vertex.x
            avg_point[1] += vertex.y
            avg_point[2] += vertex.z
        avg_point = [x / len(self.vertices) for x in avg_point]
        for vertex in self.vertices:
            vertex.x -= avg_point[0]
            vertex.y -= avg_point[1]
            vertex.z -= avg_point[2]
            norm = vertex.norm()
            vertex.x /= norm
            vertex.y /= norm
            vertex.z /= norm

    def give_mapped_mesh(self, mapping):
        vertex_list = []  # stores the vertex list of spherical mesh

        for vertex in self.vertices:
            vertex_list.append(mapping[vertex.index])

        self.vertices = vertex_list

    def get_faces_of_vertex(self, vertex):
        """get the faces of a vertex"""
        faces = []
        halfedge = vertex.halfedge.opposite.prev
        while halfedge != vertex.halfedge:
            faces.append(halfedge.facet)
            halfedge = halfedge.opposite.prev
        faces.append(halfedge.facet)
        return faces

    def gauss_map(self):
        # mapping = {}
        for vertex in self.vertices:
            face_list = self.get_faces_of_vertex(vertex)
            for face in face_list:
                face_normal = face.get_normal()
                vertex.x += face_normal[0]
                vertex.y += face_normal[1]
                vertex.z += face_normal[2]
            norm = vertex.norm()
            vertex.x /= norm
            vertex.y /= norm
            vertex.z /= norm

    def vertex_normal(self, vertex):
        face_list = self.get_faces_of_vertex(vertex)
        normal = [0.0, 0.0, 0.0]
        for face in face_list:
            normal[0] += face.get_normal()[0]
            normal[1] += face.get_normal()[1]
            normal[2] += face.get_normal()[2]
        norm = ConformalMap.norm(normal)
        normal[0] /= norm
        normal[1] /= norm
        normal[2] /= norm
        return normal

    def init_kuv(self):
        for he in self.halfedges:
            v1 = he.vertex
            v2 = he.prev.vertex
            v3 = he.next.vertex
            alpha = (
                0.5
                * ConformalMap.inner_product(
                    Vertex.subtract(v1, v3), Vertex.subtract(v2, v3)
                )
                / (
                    ConformalMap.norm(
                        ConformalMap.cross_product(
                            Vertex.subtract(v1, v3), Vertex.subtract(v2, v3)
                        )
                    )
                )
            )
            v3 = he.opposite.next.vertex
            beta = (
                0.5
                * ConformalMap.inner_product(
                    Vertex.subtract(v1, v3), Vertex.subtract(v2, v3)
                )
                / (
                    ConformalMap.norm(
                        ConformalMap.cross_product(
                            Vertex.subtract(v1, v3), Vertex.subtract(v2, v3)
                        )
                    )
                )
            )
            he.kuv = alpha + beta

    def compute_energy(self, string="tuette"):
        energy = 0.0
        explored_he = []
        for he in self.halfedges:
            if not he.opposite.index in explored_he:
                explored_he.append(he.index)
                v1 = he.vertex
                v2 = he.opposite.vertex
                v12 = Vertex.subtract(v1, v2)
                if string == "harmonic":
                    energy += he.kuv * ConformalMap.norm2(v12)
                else:
                    energy += ConformalMap.norm2(v12)

        return energy

    def get_vertices_of_vertex(self, vertex):
        """get the neigboring vertices of a vertex"""
        vertex_list = []
        halfede_list = []
        halfedge = vertex.halfedge.opposite.prev
        while halfedge != vertex.halfedge:
            halfede_list.append(halfedge)
            vertex_list.append(halfedge.opposite.vertex)
            halfedge = halfedge.opposite.prev
        vertex_list.append(halfedge.opposite.vertex)
        halfede_list.append(halfedge)
        return vertex_list, halfede_list

    def calculate_abs_dev(self, gradient, vertex):
        # point = [vertex.x, vertex.y, vertex.z]
        vertex_normal = self.vertex_normal(vertex)
        innerprod = ConformalMap.inner_product(gradient, vertex_normal)
        return [x - innerprod * y for x, y in zip(gradient, vertex_normal)]

    def compute_gradient(self, string="tuette"):
        for vertex in self.vertices:
            gradient = [0.0, 0.0, 0.0]
            neighbor_vertices, neighbor_halfedges = self.get_vertices_of_vertex(vertex)
            for i in range(len(neighbor_vertices)):
                temp_grad = Vertex.subtract(vertex, neighbor_vertices[i])
                if string == "harmonic":
                    gradient = [x + y for x, y in zip(gradient, temp_grad)]
                    gradient = [x * neighbor_halfedges[i].kuv for x in gradient]
                else:
                    gradient = [x + y for x, y in zip(gradient, temp_grad)]
            vertex.abs_gradient = self.calculate_abs_dev(gradient, vertex)

    def update_mesh(self, dt=0.05):
        for vertex in self.vertices:
            vertex.x -= dt * vertex.abs_gradient[0]
            vertex.y -= dt * vertex.abs_gradient[1]
            vertex.z -= dt * vertex.abs_gradient[2]
            norm = vertex.norm()
            vertex.x /= norm
            vertex.y /= norm
            vertex.z /= norm

    def tuette_map(self, dE=0.00001):
        curr_energy = self.compute_energy()
        print(f"Initial Tuette energy: {curr_energy}")
        prev_energy = 1000
        i = 0
        while abs(curr_energy - prev_energy) > dE:
            prev_energy = curr_energy
            self.compute_gradient()
            self.update_mesh()
            curr_energy = self.compute_energy()
            print(
                f"Current Tuette energy - {i}: {curr_energy}, Diff: {curr_energy - prev_energy}"
            )
            i += 1

    def conformal_mapping(self):
        self.init_kuv()
        self.star_map()
        # self.give_mapped_mesh(mapping)
        self.tuette_map()

    def cross_product(x1, x2):
        return [
            x1[1] * x2[2] - x1[2] * x2[1],
            x1[2] * x2[0] - x1[0] * x2[2],
            x1[0] * x2[1] - x1[1] * x2[0],
        ]

    def inner_product(x1, x2):
        return x1[0] * x2[0] + x1[1] * x2[1] + x1[2] * x2[2]

    def norm(x):
        return math.sqrt(x[0] ** 2 + x[1] ** 2 + x[2] ** 2)

    def norm2(x):
        return x[0] ** 2 + x[1] ** 2 + x[2] ** 2
