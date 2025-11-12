from pymol import cmd


def super_ply(
    pdb_query: str,
    ply_query: str,
    pdb_ref: str,
):
    """
    pdb_query: str (object name of the query pdb)
    ply_query: str (object name, the surface ply to describe pdb_query)
    pdb_ref: str (object name of the reference pdb to superimpose)
    """

    cmd.super(pdb_query, pdb_ref)
    transformation_matrix = cmd.get_object_matrix(pdb_query)
    cmd.set_object_ttt(ply_query, transformation_matrix)
