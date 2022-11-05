import halfedge_mesh


class ConformalMap(halfedge_mesh.HalfedgeMesh):
    """Conformal mapping super class of halfedge mesh"""

    def __init__(self, mesh_path):
        super().__init__(mesh_path)

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
