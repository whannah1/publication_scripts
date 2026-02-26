import matplotlib.pyplot as plt, cartopy
import matplotlib.ticker as mticker

fig_file = 'figs_sohip/sohip.locations_map.v1.png'

locations = [
    (-3.4557, -5.7297),
    (9.8802, 17.4857),
    (-49.4319, -83.461),
    (-49.185, -25.908),
    (-6.178, -56.385),
    (-52.899, 126.088),
    (-52.673, -129.993),
    (-2.744, 28.908)
]

projection = cartopy.crs.PlateCarree()

# Separate into latitude and longitude lists
latitudes = [loc[0] for loc in locations]
longitudes = [loc[1] for loc in locations]

# Create figure with map projection
fig = plt.figure(figsize=(15, 10))
ax = plt.axes(projection=projection)
ax.set_global()


# Add map features
# ax.add_feature(cartopy.feature.LAND, facecolor='lightgray')
# ax.add_feature(cartopy.feature.OCEAN, facecolor='lightblue')
ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5)
# ax.add_feature(cartopy.feature.BORDERS, linewidth=0.5, alpha=0.5)

gl = ax.gridlines(draw_labels=True, alpha=0.3, linewidth=1, color='gray')
gl.xlabels_top  = False
gl.ylabels_left = False
# gl.xlines       = False
gl.xlocator     = mticker.FixedLocator([-180, -90, 0, 90, 180])
gl.ylocator     = mticker.FixedLocator([-60, -30, 0, 30, 60])
# gl.xformatter   = LONGITUDE_FORMATTER
# gl.yformatter   = LATITUDE_FORMATTER
# gl.xlabel_style = {'size': 15, 'color': 'gray'}
# gl.xlabel_style = {'color': 'red', 'weight': 'bold'}


# Plot the locations
ax.scatter(longitudes, latitudes, c='red', s=100, marker='o', 
           edgecolors='black', linewidths=0, zorder=5,
           transform=projection)

# Add point labels
for i, (lat, lon) in enumerate(locations):
    ax.text( lon+2, lat+2, f'{i+1}', fontsize=10, transform=projection)

# ax.set_title('World Map with Location Markers', fontsize=14, fontweight='bold')


fig.savefig(fig_file, dpi=100, bbox_inches='tight')
plt.close(fig)

print(f'\n{fig_file}\n')
