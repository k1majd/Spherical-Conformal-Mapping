"""loop subdivision

author: Keyvan Majd
email: majd@asu.edu 
"""

import os
import argparse
import halfedge_mesh


def arg_parser():
    """_summary_

    Returns:
        _type_: _description_
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m",
        "--mesh",
        type=str,
        nargs="?",
        default="brain",
        help="Mesh name. Type: str.",
    )
    parser.add_argument(
        "-it",
        "--iterations",
        type=int,
        nargs="?",
        default=2,
        help="Number of iterations. Type: int.",
    )
    return parser.parse_args()


def save_halfmesh_as_obj(mesh, file_name):
    file_name = file_name + ".obj"
    with open(file_name, "w") as open_file:
        for vertex in mesh.vertices:
            lv = vertex.get_vertex()
            open_file.write("v {} {} {} \n".format(lv[0], lv[1], lv[2]))

        for face in mesh.facets:
            open_file.write("f {} {} {}\n".format(face.a + 1, face.b + 1, face.c + 1))


def save_halfmesh_as_off(mesh, file_name):
    file_name = file_name + ".off"
    with open(file_name, "w") as open_file:
        open_file.write("OFF\n")
        open_file.write(f"{len(mesh.vertices)} {len(mesh.facets)} 0\n")
        for vertex in mesh.vertices:
            lv = vertex.get_vertex()
            open_file.write("{} {} {} \n".format(lv[0], lv[1], lv[2]))

        for face in mesh.facets:
            open_file.write("3 {} {} {}\n".format(face.a, face.b, face.c))


def main(mesh_name):
    cwd = os.path.dirname(os.path.realpath(__file__))
    # mesh_path = f"halfedge_mesh/tests/data/{mesh_name}.off"
    mesh_path = cwd + f"/tests/data/{mesh_name}.off"
    save_path = f"{mesh_name}1"

    mesh = halfedge_mesh.ConformalMap(mesh_path)
    mesh.conformal_mapping()
    save_halfmesh_as_obj(mesh, f"{save_path}")

    # for i in range(iterations - 1):
    #     print(f"loop subdivision, iteration {i + 1}")
    #     mesh = halfedge_mesh.ConformalMap(mesh_path)
    #     if i > 0:
    #         # remove the old mesh
    #         os.remove(mesh_path)
    #     new_mesh = mesh.conformal_mapping()
    #     save_halfmesh_as_off(new_mesh, f"{save_path}_{i+1}")
    #     mesh_path = f"{save_path}_{i+1}.off"
    # print("loop subdivision, iteration", iterations)
    # mesh = halfedge_mesh.SubDivider(mesh_path)
    # new_mesh = mesh.loop_subdivision()
    # if iterations > 1:
    #     os.remove(mesh_path)
    # save_halfmesh_as_obj(new_mesh, f"{save_path}_{iterations}")
    # current_directory = os.getcwd() + "/" + f"{save_path}_{iterations}.obj"
    # print("mesh saved at", current_directory)


if __name__ == "__main__":
    args = arg_parser()
    main(args.mesh)
