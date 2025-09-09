## Ground track plot from orbital parameters


import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.animation import FuncAnimation

# ============================================================
# Constants
# ============================================================
pi = np.pi
mu = 3.986004418e14  # [m^3/s^2] Earth GM
rEarth = 6371000     # mean Earth radius [m]

# WGS84 ellipsoid parameters
a = 6378137.0        # semi-major axis [m]
f = 1 / 298.257223563
e2 = f * (2 - f)     # eccentricity squared

# ====================================================================================
# Rotation Matrices: M313 for ECI to otìrbital reference frame, and M3 for ECI to ECEF
# ====================================================================================
def M313(raan, inc, aop):
    # 313 rotation (ECI -> orbital plane)
    M3raan = np.array([
        [np.cos(raan), -np.sin(raan), 0],
        [np.sin(raan),  np.cos(raan), 0],
        [0, 0, 1]
    ])
    M1inc = np.array([
        [1, 0, 0],
        [0, np.cos(inc), -np.sin(inc)],
        [0, np.sin(inc),  np.cos(inc)]
    ])
    M3aop = np.array([
        [np.cos(aop), -np.sin(aop), 0],
        [np.sin(aop),  np.cos(aop), 0],
        [0, 0, 1]
    ])
    return M3aop @ M1inc @ M3raan

def M3alphaGr(alphaGr):
    #rotation around Z axis (Earth rotation axis) ECI -> ECEF
    return np.array([
        [np.cos(alphaGr), -np.sin(alphaGr), 0],
        [np.sin(alphaGr),  np.cos(alphaGr), 0],
        [0, 0, 1]
    ])

# ============================================================
# Coordinate Conversions
# ============================================================
def ECEFtoGeodetic(x, y, z):
    #Conversion from ECEF to lat/long Geodetic coordinates (WGS84)
    lon = np.arctan2(y, x)
    r = np.sqrt(x**2 + y**2)
    lat = np.arctan2(z, r * (1 - e2))  # approximate value, an iterative method would be more accurate
    return lon, lat

# ============================================================
# Ground Track Computation
# ============================================================
def compute_ground_track(sma, ecc, inc, raan, aop, n_points=360, start_time=0):
    # The actual functions that converts orbital parameters into the ground track
    # Orbital parameters
    p = sma * (1 - ecc**2)
    ta = np.linspace(0, 2*pi, n_points)
    r_orb = p / (1 + ecc*np.cos(ta))

    # Position vectors in orbital plane
    r_orb_vec = np.vstack((r_orb*np.cos(ta),
                           r_orb*np.sin(ta),
                           np.zeros_like(ta)))

    # Transform to ECI
    Morbtoeci = np.linalg.inv(M313(raan, inc, aop))
    r_eci = Morbtoeci @ r_orb_vec

    # Time evolution
    T_orbit = np.sqrt(sma**3 / mu)
    omega_earth = 7.2921159e-5  # rad/s
    dt = T_orbit / n_points
    times = np.arange(n_points) * dt
    alphaGr = (start_time/24)*2*pi + omega_earth * times

    # Transform to ECEF
    r_ecef = np.zeros_like(r_eci)
    for k in range(n_points):
        r_ecef[:, k] = M3alphaGr(alphaGr[k]) @ r_eci[:, k]

    # Convert to geodetic coordinates
    lon, lat = np.zeros(n_points), np.zeros(n_points)
    for k in range(n_points):
        lon[k], lat[k] = ECEFtoGeodetic(*r_ecef[:, k])

    return np.degrees(lat), np.degrees(lon)

# ============================================================
# Plotting
# ============================================================
def plot_ground_track(lat, lon):
    #Static plot of the ground track
    fig = plt.figure(figsize=(10, 5))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_global()
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.gridlines(draw_labels=True)
    ax.plot(lon, lat, 'r-', transform=ccrs.Geodetic())
    plt.title("Satellite Ground Track")
    plt.show()

def animate_ground_track(lat, lon):
    #Animated ground track with moving satellite
    fig = plt.figure(figsize=(10, 5))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_global()
    ax.coastlines()
    ax.gridlines()

    line, = ax.plot([], [], 'r-', transform=ccrs.Geodetic())
    point, = ax.plot([], [], 'bo', transform=ccrs.Geodetic())

    def init():
        line.set_data([], [])
        point.set_data([], [])
        return line, point

    def update(frame):
        line.set_data(lon[:frame], lat[:frame])
        point.set_data(lon[frame], lat[frame])
        return line, point

    ani = FuncAnimation(fig, update, frames=len(lat), init_func=init,
                        interval=50, blit=True)
    plt.show()

# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    # Example orbital parameters
    sma = (6378 + 600) * 1000  # semi-major axis [m]
    ecc = 0.01
    inc = np.radians(60)
    raan = np.radians(20)
    aop = 0

    # Compute ground track
    lat, lon = compute_ground_track(sma, ecc, inc, raan, aop,
                                    n_points=720, start_time=15.5)

    # Plot results
    plot_ground_track(lat, lon)
    animate_ground_track(lat, lon)

