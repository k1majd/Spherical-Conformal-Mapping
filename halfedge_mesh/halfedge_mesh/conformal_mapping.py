import halfedge_mesh


class ConformalMap(halfedge_mesh.HalfedgeMesh):
    """Conformal mapping super class of halfedge mesh"""

    def __init__(self, mesh_path):
        super().__init__(mesh_path)
        self.boundary_vertex_set = self._collect_boundary_vertices()

    def if_boundary_halfedge(self, halfedge):
        """checks if the halfedge is a boundary halfedge"""
        if halfedge.opposite is None:
            return True
        else:
            return False

    def if_boundary_vertex(self, vertex):
        """checks if the vertex is a boundary vertex"""
        if vertex in self.boundary_vertex_set:
            return True
        else:
            return False

    def even_vertex_adjacency(self, vertex):
        """find the adjacent vertices of an even vertex"""
        vertex_neighbor_list = [vertex]
        if self.if_boundary_vertex(vertex):  # boundary vertex
            halfedge = vertex.halfedge
            while halfedge.opposite is not None:
                halfedge = halfedge.opposite.prev
            vertex_neighbor_list.append(halfedge.prev.vertex)
            halfedge = vertex.halfedge
            while halfedge.next.opposite is not None:
                halfedge = halfedge.next.opposite
            vertex_neighbor_list.append(halfedge.next.vertex)
        else:  # interior vertex
            halfedge = vertex.halfedge.opposite.prev
            while halfedge != vertex.halfedge:
                vertex_neighbor_list.append(halfedge.opposite.vertex)
                halfedge = halfedge.opposite.prev
            vertex_neighbor_list.append(halfedge.opposite.vertex)
        return vertex_neighbor_list

    def odd_vertex_adjacency(self, halfedge):
        """find the adjacent vertices of an odd vertex"""
        vertex_neighbor_list = []

        if self.if_boundary_halfedge(halfedge):  # boundary halfedge
            vertex_neighbor_list.append(halfedge.vertex)
            vertex_neighbor_list.append(halfedge.next.next.vertex)
        else:  # interior halfedge
            vertex_neighbor_list.append(halfedge.vertex)
            vertex_neighbor_list.append(halfedge.opposite.vertex)
            vertex_neighbor_list.append(halfedge.next.vertex)
            vertex_neighbor_list.append(halfedge.opposite.next.vertex)

        return vertex_neighbor_list

    def compute_even_vertex(self, vertex, index=None):
        """compute the even vertex"""
        vertex_neighbor_list = self.even_vertex_adjacency(vertex)
        n = len(vertex_neighbor_list)
        even_vertex = halfedge_mesh.Vertex(index=index)
        if n == 3:  # boundary vertex
            even_vertex.x = (
                6 * vertex_neighbor_list[0].x
                + vertex_neighbor_list[1].x
                + vertex_neighbor_list[2].x
            ) / 8
            even_vertex.y = (
                6 * vertex_neighbor_list[0].y
                + vertex_neighbor_list[1].y
                + vertex_neighbor_list[2].y
            ) / 8
            even_vertex.z = (
                6 * vertex_neighbor_list[0].z
                + vertex_neighbor_list[1].z
                + vertex_neighbor_list[2].z
            ) / 8
        else:  # interior vertex
            k = n - 1
            beta = 3 / (8 * n) if n != 3 else 3 / 16
            even_vertex.x = (1 - k * beta) * vertex_neighbor_list[0].x + beta * sum(
                [vertex.x for vertex in vertex_neighbor_list[1:]]
            )
            even_vertex.y = (1 - k * beta) * vertex_neighbor_list[0].y + beta * sum(
                [vertex.y for vertex in vertex_neighbor_list[1:]]
            )
            even_vertex.z = (1 - k * beta) * vertex_neighbor_list[0].z + beta * sum(
                [vertex.z for vertex in vertex_neighbor_list[1:]]
            )
        return even_vertex

    def compute_odd_vertex(self, halfedge, index=None):
        """compute the odd vertex"""
        vertex_neighbor_list = self.odd_vertex_adjacency(halfedge)
        n = len(vertex_neighbor_list)
        odd_vertex = halfedge_mesh.Vertex(index=index)
        if n == 2:  # boundary vertex
            odd_vertex.x = (vertex_neighbor_list[0].x + vertex_neighbor_list[1].x) / 2
            odd_vertex.y = (vertex_neighbor_list[0].y + vertex_neighbor_list[1].y) / 2
            odd_vertex.z = (vertex_neighbor_list[0].z + vertex_neighbor_list[1].z) / 2
        elif n == 4:  # interior vertex
            odd_vertex.x = (
                (vertex_neighbor_list[0].x + vertex_neighbor_list[1].x) * 3
                + (vertex_neighbor_list[2].x + vertex_neighbor_list[3].x)
            ) / 8
            odd_vertex.y = (
                (vertex_neighbor_list[0].y + vertex_neighbor_list[1].y) * 3
                + (vertex_neighbor_list[2].y + vertex_neighbor_list[3].y)
            ) / 8
            odd_vertex.z = (
                (vertex_neighbor_list[0].z + vertex_neighbor_list[1].z) * 3
                + (vertex_neighbor_list[2].z + vertex_neighbor_list[3].z)
            ) / 8
        else:
            raise ValueError("odd vertex should have 2 or 4 neighbors")
        return odd_vertex

    def loop_subdivision(self):
        """perform loop subdivision on the mesh"""
        mesh_new = halfedge_mesh.HalfedgeMesh()
        vertex_list = []  # stores the vertex list of new mesh
        face_list = []  # stores the face list of new mesh
        face_index = 0  # index of current created face
        vertex_index = 0  # index of current created vertex
        odd_explored = {}  # map each explored halfedge to its odd vertex
        even_explored = {}  # map each explored vertex to its even vertex
        vertex_index_recorder = {
            0: None,
            1: None,
            2: None,
            3: None,
            4: None,
            5: None,
        }  # record the vertices of a face
        for face in self.facets:
            vertex_tracker = 0  # track the verteces of a face
            halfedge = face.halfedge
            for _ in range(3):
                # each half edge corresponds to a pair of
                # 1. odd and 2. even vertex

                # 1. create odd vertex
                if self.if_boundary_halfedge(halfedge):
                    odd_vertex = self.compute_odd_vertex(halfedge, index=vertex_index)
                    vertex_list.append(odd_vertex)
                    odd_explored[halfedge.index] = odd_vertex
                    vertex_index_recorder[vertex_tracker] = odd_vertex.index
                    vertex_index += 1
                    vertex_tracker += 1
                else:
                    if (
                        halfedge.opposite.index not in odd_explored.keys()
                    ):  # if the odd vertex is not created
                        odd_vertex = self.compute_odd_vertex(
                            halfedge, index=vertex_index
                        )
                        vertex_list.append(odd_vertex)
                        odd_explored[halfedge.index] = odd_vertex
                        vertex_index_recorder[vertex_tracker] = odd_vertex.index
                        vertex_index += 1
                        vertex_tracker += 1
                    else:  # if the odd vertex is created
                        odd_vertex = odd_explored[halfedge.opposite.index]
                        vertex_index_recorder[vertex_tracker] = odd_vertex.index
                        vertex_tracker += 1

                # 2. create even vertex
                if (
                    halfedge.vertex.index not in even_explored.keys()
                ):  # if the even vertex is not created
                    even_vertex = self.compute_even_vertex(
                        halfedge.vertex, index=vertex_index
                    )
                    vertex_list.append(even_vertex)
                    even_explored[halfedge.vertex.index] = even_vertex
                    vertex_index_recorder[vertex_tracker] = even_vertex.index
                    vertex_index += 1
                    vertex_tracker += 1
                else:  # if the even vertex is created
                    even_vertex = even_explored[halfedge.vertex.index]
                    vertex_index_recorder[vertex_tracker] = even_vertex.index
                    vertex_tracker += 1

                halfedge = halfedge.next

            # create the new faces
            face_list.append(
                halfedge_mesh.Facet(
                    a=vertex_index_recorder[0],
                    b=vertex_index_recorder[1],
                    c=vertex_index_recorder[2],
                    index=face_index,
                )
            )
            face_index += 1
            face_list.append(
                halfedge_mesh.Facet(
                    a=vertex_index_recorder[2],
                    b=vertex_index_recorder[3],
                    c=vertex_index_recorder[4],
                    index=face_index,
                )
            )
            face_index += 1
            face_list.append(
                halfedge_mesh.Facet(
                    a=vertex_index_recorder[4],
                    b=vertex_index_recorder[5],
                    c=vertex_index_recorder[0],
                    index=face_index,
                )
            )
            face_index += 1
            face_list.append(
                halfedge_mesh.Facet(
                    a=vertex_index_recorder[0],
                    b=vertex_index_recorder[2],
                    c=vertex_index_recorder[4],
                    index=face_index,
                )
            )
            face_index += 1

        mesh_new.vertices = vertex_list
        mesh_new.facets = face_list
        return mesh_new

    def star_map(self):
        """map the vertices in the star of a vertex"""
        mapping = {}
        avg_point = [0.0, 0.0, 0.0]
        for vertex in self.vertices:
            avg_point[0] += vertex.x
            avg_point[1] += vertex.y
            avg_point[2] += vertex.z
        avg_point = [x / len(self.vertices) for x in avg_point]
        for vertex in self.vertices:
            v_temp = halfedge_mesh.Vertex(
                index=vertex.index,
                x=vertex.x,
                y=vertex.y,
                z=vertex.z,
                halfedge=vertex.halfedge,
            )
            v_temp.x -= avg_point[0]
            v_temp.y -= avg_point[1]
            v_temp.z -= avg_point[2]
            norm = v_temp.norm()
            v_temp.x /= norm
            v_temp.y /= norm
            v_temp.z /= norm
            mapping[vertex.index] = v_temp

        return mapping

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
        mapping = {}
        for vertex in self.vertices:
            new_vertex = halfedge_mesh.Vertex(
                index=vertex.index, halfedge=vertex.halfedge
            )
            face_list = self.get_faces_of_vertex(vertex)
            for face in face_list:
                face_normal = face.get_normal()
                new_vertex.x += face_normal[0]
                new_vertex.y += face_normal[1]
                new_vertex.z += face_normal[2]
            norm = new_vertex.norm()
            new_vertex.x /= norm
            new_vertex.y /= norm
            new_vertex.z /= norm
            mapping[vertex.index] = new_vertex
        return mapping

    def conformal_mapping(self):
        mapping = self.star_map()
        self.give_mapped_mesh(mapping)
        print("here")

    def _collect_boundary_vertices(self):
        """find the boundary vertices of the mesh"""
        boundary_halfedges = []
        for halfedge in self.halfedges:
            if self.if_boundary_halfedge(halfedge):
                boundary_halfedges.append(halfedge.vertex)
        return set(boundary_halfedges)
