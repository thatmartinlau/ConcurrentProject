import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Load telemetry data
df = pd.read_csv("telemetry.csv")

# Get number of bodies and steps
num_bodies = df['body'].nunique()
num_steps = df['step'].nunique()

# Assign colors if needed
colors = ['yellow', 'blue', 'red', 'green', 'purple', 'orange']

fig, ax = plt.subplots()
scat = ax.scatter([], [], s=100)

def init():
    ax.set_xlim(df['x'].min() - 10, df['x'].max() + 10)
    ax.set_ylim(df['y'].min() - 10, df['y'].max() + 10)
    return scat,

def update(frame):
    frame_data = df[df['step'] == frame]
    scat.set_offsets(frame_data[['x', 'y']])
    scat.set_color([colors[i % len(colors)] for i in frame_data['body']])
    ax.set_title(f"Step {frame}")
    return scat,

ani = animation.FuncAnimation(fig, update, frames=num_steps,
                              init_func=init, blit=True, interval=50)

ani.save("telemetry_animation.gif", writer="pillow")
plt.show()
