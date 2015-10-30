from fidimag.common import CuboidMesh


def test_mesh1():
    mesh = CuboidMesh(nx=5, ny=3, nz=2, dx=0.23, dy=0.41)
    assert len(mesh.coordinates) == 5 * 3 * 2
    assert mesh.ny * mesh.nz == 6
    assert tuple(mesh.coordinates[mesh.index(0, 0, 0)]) == (0.23 / 2, 0.41 / 2, 0.5)
    assert tuple(mesh.coordinates[mesh.index(3, 2, 1)]) == (
        (3 + 0.5) * 0.23, (2 + 0.5) * 0.41, 1 + 0.5)

if __name__ == '__main__':
    test_mesh1()
