import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

# Read the CSV file
data = pd.read_csv('telemetry.csv')

# Create two visualizations: static plot and animation
plt.style.use('dark_background')  # Optional: use dark theme for space-like appearance

# 1. Static Plot (shows complete trajectories)
plt.figure(figsize=(12, 12))

# Plot trajectories
plt.plot(data['Sun_x'], data['Sun_y'], 'yo-', label='Sun', markersize=10)
plt.plot(data['Mercury_x'], data['Mercury_y'], 'mo-', label='Mars', markersize=6)
plt.plot(data['Venus_x'], data['Venus_y'], 'go-', label='Mars', markersize=6)
plt.plot(data['Earth_x'], data['Earth_y'], 'bo-', label='Earth', markersize=6)
plt.plot(data['Mars_x'], data['Mars_y'], 'ro-', label='Mars', markersize=6)


plt.title('Solar System Simulation Trajectories')
plt.xlabel('X Position (m)')
plt.ylabel('Y Position (m)')
plt.grid(True, alpha=0.3)
plt.axis('equal')  # Make sure orbits appear circular
plt.legend()

# Scientific notation for axis labels
plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))

# Save static plot
plt.savefig('trajectories.png')
plt.show()

# 2. Animated Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))
fig.suptitle('Solar System Simulation')

# First subplot: Complete view
def animate(frame):
    ax1.clear()
    ax2.clear()
    
    # Complete view (first subplot)
    ax1.plot(data['Sun_x'][:frame], data['Sun_y'][:frame], 'yo-', label='Sun', markersize=10)
    ax1.plot(data['Mercury_x'][:frame], data['Mercury_y'][:frame], 'mo-', label='Mercury', markersize=10)
    ax1.plot(data['Venus_x'][:frame], data['Venus_y'][:frame], 'go-', label='Venus', markersize=10)
    ax1.plot(data['Earth_x'][:frame], data['Earth_y'][:frame], 'bo-', label='Earth', markersize=6)
    ax1.plot(data['Mars_x'][:frame], data['Mars_y'][:frame], 'ro-', label='Mars', markersize=6)
    
    
    ax1.set_title(f'Time: {data["time"][frame]:.1f} seconds')
    ax1.grid(True, alpha=0.3)
    ax1.axis('equal')
    ax1.legend()
    
    # Sun-centered zoomed view (second subplot)
    zoom_range = 2e11  # adjust based on your data
    ax2.plot(data['Sun_x'][frame], data['Sun_y'][frame], 'yo', label='Sun', markersize=10)
    ax2.plot(data['Mercury_x'][frame], data['Mercury_y'][frame], 'mo', label='Mercury', markersize=10)
    ax2.plot(data['Venus_x'][frame], data['Venus_y'][frame], 'go', label='Venus', markersize=10)
    ax2.plot(data['Earth_x'][frame], data['Earth_y'][frame], 'bo', label='Earth', markersize=6)
    ax2.plot(data['Mars_x'][frame], data['Mars_y'][frame], 'ro', label='Mars', markersize=6)
    
    ax2.set_xlim(-zoom_range, zoom_range)
    ax2.set_ylim(-zoom_range, zoom_range)
    ax2.set_title('Zoomed View')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # Scientific notation for both plots
    ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0))

anim = FuncAnimation(fig, animate, frames=len(data), interval=50)

# Save animation
anim.save('solar_system.gif', writer='pillow')
plt.show()

# Print some statistics
bodies = ['Sun', 'Mercury', 'Venus','Earth', 'Mars']
for body in bodies:
    print(f"\n{body} statistics:")
    print(f"X range: {data[f'{body}_x'].min():.2e} to {data[f'{body}_x'].max():.2e}")
    print(f"Y range: {data[f'{body}_y'].min():.2e} to {data[f'{body}_y'].max():.2e}")
    
# Calculate and print orbital velocities
for body in ['Mercury', 'Venus', 'Earth', 'Mars']:
    dx = np.diff(data[f'{body}_x'])
    dy = np.diff(data[f'{body}_y'])
    dt = np.diff(data['time'])
    velocities = np.sqrt(dx**2 + dy**2) / dt
    avg_velocity = np.mean(velocities)
    print(f"\n{body} average orbital velocity: {avg_velocity:.2e} m/s")