import numpy as np
# Monoclinic cell parameters
a = 8.5550
b = 5.9680
c = 20.0331
beta_deg = 110.055
# Convert beta to radians
beta_rad = np.deg2rad(beta_deg)
sin_beta = np.sin(beta_rad)
cos_beta = np.cos(beta_rad)
# Miller indices of the given plane
h1, k1, l1 = 1, 0, 4
def scalar_product(h1, k1, l1, h2, k2, l2, a, b, c, sin_beta, cos_beta):
    term1 = (h1 * h2) / (a**2 * sin_beta**2)
    term2 = (k1 * k2) / (b**2)
    term3 = (l1 * l2) / (c**2 * sin_beta**2)
    term4 = cos_beta / (a * c * sin_beta**2) * (h1 * l2 + h2 * l1)
    return term1 + term2 + term3 - term4
# Calculate coefficients for the linear relationship
A = 1/(a**2 * sin_beta**2) - 4*cos_beta/(a*c*sin_beta**2)
B = 4/(c**2 * sin_beta**2) - cos_beta/(a*c*sin_beta**2)
best_S = None
best_indices = None
# Search for the best integer solution for l2 in a reasonable range
for l2 in range(1, 30):
    h2 = -B/A * l2
    h2_int = int(round(h2))
    S = scalar_product(h1, k1, l1, h2_int, 0, l2, a, b, c, sin_beta, cos_beta)
    if best_S is None or abs(S) < abs(best_S):
        best_S = S
        best_indices = (h2_int, 0, l2)
    print(f"Trying (h2, k2, l2) = ({h2_int}, 0, {l2}): S = {S}")
print(f"\nBest indices: {best_indices}, S = {best_S}")