import numpy as np
import os

# Create directory structure
os.makedirs("data/mapping/panda/100", exist_ok=True)
os.makedirs("data/mapping/panda/10_new", exist_ok=True)

print("Creating corrected bootstrap data...")

# Generate 100 random 3x3 positive definite matrices (teacher lookup data)
np.random.seed(42)
teacher_data = []
for i in range(100):
    A = np.random.randn(3, 3)
    A = A @ A.T + np.eye(3) * 0.1  # Ensure positive definite
    teacher_data.append(A.flatten())

teacher_data = np.array(teacher_data)
np.savetxt("data/mapping/panda/100/manipulabilities.csv", teacher_data, delimiter=",")

# Generate 20 interpolated points WITH TIME COLUMN (10 columns total)
np.random.seed(123)
interpolated_data = []
for i in range(20):
    time_step = i * 0.1  # Time column
    A = np.random.randn(3, 3)
    A = A @ A.T + np.eye(3) * 0.2  # Different seed for variety
    row = np.concatenate([[time_step], A.flatten()])  # 1 + 9 = 10 columns
    interpolated_data.append(row)

interpolated_data = np.array(interpolated_data)
np.savetxt("data/mapping/panda/10_new/manipulabilities_interpolated.csv", interpolated_data, delimiter=",")

print("Files created:")
print("  data/mapping/panda/100/manipulabilities.csv (100x9)")
print("  data/mapping/panda/10_new/manipulabilities_interpolated.csv (20x10)")
print("\nNow you can run:")
print("  python3 py/generate_artificial_data_r_s_t.py data/mapping 100 10_new")

