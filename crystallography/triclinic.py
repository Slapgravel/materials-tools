import numpy as np

def cartesian_to_miller(a, b, c, r):
    """
    Convert a direction vector in Cartesian coordinates to Miller indices.

    Parameters:
    - a (list): Lattice vector a (in Cartesian coordinates)
    - b (list): Lattice vector b (in Cartesian coordinates)
    - c (list): Lattice vector c (in Cartesian coordinates)
    - r (list): Direction vector in Cartesian coordinates

    Returns:
    - miller_indices (list): Miller indices [u, v, w] as integers
    """
    # Create lattice vectors matrix
    A = np.array([a, b, c]).T

    # Calculate inverse of A
    A_inv = np.linalg.inv(A)

    # Solve for u, v, w
    uvw = np.dot(A_inv, r)

    # Round to nearest integers
    miller_indices = [round(x) for x in uvw]

    return miller_indices

# Example usage:
lattice_vectors = {
    'a': [4.23, 0, 0],
    'b': [-2.12, 6.78, 0],
    'c': [0, -0.34, 3.45]
}

direction_vector = [1.23, 4.56, 7.89]

miller_indices = cartesian_to_miller(**lattice_vectors, r=direction_vector)
print("Miller Indices:", miller_indices)