import numpy as np
import matplotlib.pyplot as plt

# Load data from CSV file
data = np.loadtxt('res.csv', delimiter=',', skiprows=1)
L = data[:, 0]  # Lattice constants
P = data[:, 1]  # Potential energy values

# Fit a second-degree polynomial to the data
coeffs = np.polyfit(L, P, 2)
poly = np.poly1d(coeffs)

# Find the optimal L by finding the minimum of the fitted polynomial
L_opt = -coeffs[1] / (2 * coeffs[0])
P_min = poly(L_opt)

# Plot the data and the fitted curve
plt.scatter(L, P, label='Data points', color='blue')
L_fit = np.linspace(min(L), max(L), 100)
P_fit = poly(L_fit)
plt.plot(L_fit, P_fit, label=f'2nd degree fit', color='red')

# Mark the optimal point on the plot
plt.plot(L_opt, P_min, 'go', label=f'Optimal L={L_opt:.3f}, P_min={P_min:.3f}')

# Labels and title
plt.xlabel('Lattice Constant (L)')
plt.ylabel('Potential Energy (P)')
plt.title('Lattice Constant vs. Potential Energy')
plt.legend()
plt.grid(True)
plt.savefig('fit.png')
# plt.show()

# Print results
print(f"Optimal Lattice Constant (L_opt): {L_opt:.3f}")
print(f"Minimum Potential Energy (P_min): {P_min:.3f}")
