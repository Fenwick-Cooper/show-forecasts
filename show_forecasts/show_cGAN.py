# Load and plot cGAN forecast data.

# To do:
#   Change "mm h**-1" to "mm/h" in the data.
#   Store and extract the model name from the data.
#   Add the initialisation time to the title.

import numpy as np
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
from cartopy.feature import ShapelyFeature
import matplotlib.pyplot as plt
from matplotlib import colors  # For consistency with Harris et. al 2022
from datetime import timedelta
import xarray as xr
from .data_utils import *


# Load a 24 hour mean forecast at a lead time of 30 to 54 hours
# To be clear about the dates use print_forecast_info(forecast_init_date)
# Arguments
#    forecast_init_date   - A datetime.datetime corresponding to when the forecast was initialised.
#    data_dir             - Directory where the data is stored.
# Returns
#    An xarray DataSet containing the cGAN rainfall forecasts.
def load_GAN_forecast(forecast_init_date, data_dir):
    d = forecast_init_date  # Shorthand
    file_name = f"{data_dir}/GAN_{d.year}{d.month:02}{d.day:02}.nc"
    data = xr.open_dataset(file_name)
    return data


# Sorts the data along the ensemble axis. Used in percentile plots.
# Arguments
#   data                    - An xarray DataSet containing the cGAN rainfall forecasts.
# Returns
#   data_sorted             - A xarray DataSet containing the cGAN rainfall forecasts sorted along the
#                             ensemble axis.
def sort_along_ensemble_axis(data):
    
    # Create a copy of the Dataset to hold the sorted data
    data_sorted = data.copy()

    # Find the axis to sort along
    member_axis = data['precipitation'].get_axis_num('member')

    # Sort values along the chosen axis
    precip_sorted = np.sort(data["precipitation"], axis=member_axis)

    # Copy the sorted data into the data_sorted Dataset
    data_sorted["precipitation"].values = precip_sorted

    return data_sorted


# Plot the ensemble mean and ensemble standard deviation of the cGAN forecast data
# at each valid time.
# Arguments
#   data                        - An xarray DataSet containing the cGAN rainfall forecasts.
#   accumulation_time='6h'      - Can be '6h', or '24h'
#   valid_time_start_hour='all' - The hour the valid time starts at. Can either be 6, 12, 18
#                                 or 0 UTC, or specify 'all' to make all four plots.
#   style=None                  - Options: 'ICPAC', 'ICPAC_heavy', 'KMD', 'EMI', 'EMI_heavy'
#   plot_units='mm/h'           - Can be 'mm/h' (default), 'mm/6h', 'mm/day' or 'mm/week'
#   region='ICPAC'              - can be 'ICPAC', 'Kenya', 'South Sudan', 'Rwanda', 'Burundi',
#                                 'Djibouti', 'Eritrea', 'Ethiopia', 'Sudan', 'Somalia',
#                                 'Tanzania', 'Uganda'
#   file_name=None              - If a file name, ending in '.png', '.jpg' or '.pdf' is
#                                 specified, the plot is saved in that format. If
#                                 valid_time_start_hour = 'all', the hour is appended to the
#                                 file name.
def plot_GAN_forecast(data, accumulation_time='6h', valid_time_start_hour='all',
                      style=None, plot_units='mm/h', region='ICPAC', file_name=None):

    # Get the units to use for plotting
    plot_norm, plot_units = get_plot_normalisation(plot_units)
    
    # To be consistent with the Harris et. al paper.
    value_range_precip = (0.1, 15 * plot_norm)
    
    # Use a style other than the default
    if (style is not None):
        plot_levels, plot_colours = get_contour_levels(style)
    
    # Load the border shapefile
    reader = shpreader.Reader("show_forecasts/shapes/GHA_shapes/gha.shp")
    borders_feature = ShapelyFeature(reader.geometries(), ccrs.PlateCarree(), facecolor='none')
    
    if (region != 'ICPAC'):
        
        # Load the regions shapefile
        reader = shpreader.Reader(f"show_forecasts/shapes/{region}_shapes/{region}_region.cpg")
        regions_feature = ShapelyFeature(reader.geometries(), ccrs.PlateCarree(), facecolor='none')
        
        # Get the extent of the region that we are looking at
        region_extent = get_region_extent(region, border_size=0.5)

    if (accumulation_time == '6h'):
        
        # Change valid_time_start_hour into the valid_time_idx
        if (str(valid_time_start_hour) == '6'):
            valid_time_idx_list = [0]
        elif (str(valid_time_start_hour) == '12'):
            valid_time_idx_list = [1]
        elif (str(valid_time_start_hour) == '18'):
            valid_time_idx_list = [2]
        elif (str(valid_time_start_hour) == '0'):
            valid_time_idx_list = [3]
        elif (valid_time_start_hour == 'all'):
            valid_time_idx_list = np.arange(len(data['valid_time']))
        else:
            print("ERROR: valid_time_start_hour must be 6, 12, 18 or 0 hours UTC or 'all'.")
            return

        valid_time_delta = timedelta(hours=6)

    elif (accumulation_time == '24h'):

        if (str(valid_time_start_hour) != '0') and (valid_time_start_hour != 'all'):
            print("ERROR: valid_time_start_hour must be 0 when accumulation_time is '24h'.")
            return
        
        valid_time_idx_list = [0]
        valid_time_delta = timedelta(hours=24)

    else:
        print("ERROR: accumulation_time must be either '6h' or '24h'.")
        return

    # There are plots for each valid time
    for valid_time_idx in valid_time_idx_list:

        # Convert the forecast valid time to a datetime.datetime format
        valid_time = datetime64_to_datetime(data['fcst_valid_time'][0,valid_time_idx].values)

        # Define the figure and each axes for the rows and columns
        fig, axs = plt.subplots(nrows=1, ncols=2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(10,5))

        # axs is a 2 dimensional array of `GeoAxes`. Flatten it into a 1-D array
        axs=axs.flatten()

        ax=axs[0]  # First plot (left)
        ax.add_feature(cfeature.COASTLINE, linewidth=1)  # Draw some features to see where we are
        if (region != 'ICPAC'):
            if (region != 'Uganda'):  # Uganda counties are too complicated
                ax.add_feature(regions_feature, linestyle=':')
            ax.set_extent(region_extent, crs=ccrs.PlateCarree())
        ax.add_feature(borders_feature)  # The borders
        ax.add_feature(cfeature.LAKES, linewidth=1,linestyle='-',edgecolor='dimgrey',facecolor='none')

        # Either plot 6h data or 24h data
        if (accumulation_time == '6h'):
            data_to_plot = np.mean(data['precipitation'][0,:,valid_time_idx,:,:], axis=0) * plot_norm
        elif (accumulation_time == '24h'):
            data_to_plot = np.mean(data['precipitation'][0,:,:,:,:], axis=(0,1)) * plot_norm

        # Actually make the plot
        if (style == None):
            c = ax.pcolormesh(data['longitude'], data['latitude'], data_to_plot,
                              norm=colors.LogNorm(*value_range_precip), cmap='YlGnBu', transform=ccrs.PlateCarree())
            cb = plt.colorbar(c, fraction=0.04)  # Add a colorbar with a nice size
        else:
            c = ax.contourf(data['longitude'], data['latitude'], data_to_plot,
                            colors=plot_colours, levels=plot_levels*plot_norm, transform=ccrs.PlateCarree())
            cb = plt.colorbar(c, fraction=0.04)  # Add a colorbar with a nice size
            cb_labels = np.round(plot_levels*plot_norm,1).astype(str).tolist()
            cb_labels[-1] = ''  # Remove the final label
            cb.set_ticks(ticks=plot_levels*plot_norm, labels=cb_labels)
        cb.set_label(f'Rainfall ({plot_units})')  # Label the colorbar
        ax.set_title(f"Ensemble mean",size=14)  # This plot's title

        # Either plot 6h data or 24h data
        if (accumulation_time == '6h'):
            data_to_plot = np.std(data['precipitation'][0,:,valid_time_idx,:,:], axis=0, ddof=1) * plot_norm
        elif (accumulation_time == '24h'):
            data_to_plot = np.sqrt(np.mean(np.var(data['precipitation'][0,:,:,:,:], axis=0, ddof=1), axis=0)) * plot_norm

        ax=axs[1]  # Second plot (right)
        ax.add_feature(cfeature.COASTLINE, linewidth=1)  # Draw some features to see where we are
        if (region != 'ICPAC'):
            if (region != 'Uganda'):  # Uganda counties are too complicated
                ax.add_feature(regions_feature, linestyle=':')
            ax.set_extent(region_extent, crs=ccrs.PlateCarree())
        ax.add_feature(borders_feature)  # The borders
        ax.add_feature(cfeature.LAKES, linewidth=1,linestyle='-',edgecolor='dimgrey',facecolor='none')
        # Actually make the plot
        if (style == None):
            c = ax.pcolormesh(data['longitude'], data['latitude'], data_to_plot,
                              norm=colors.LogNorm(*value_range_precip), cmap='YlGnBu', transform=ccrs.PlateCarree())
            cb = plt.colorbar(c, fraction=0.04)  # Add a colorbar with a nice size
        else:
            c = ax.contourf(data['longitude'], data['latitude'], data_to_plot,
                            colors=plot_colours, levels=plot_levels*plot_norm, transform=ccrs.PlateCarree())
            cb = plt.colorbar(c, fraction=0.04)  # Add a colorbar with a nice size
            cb_labels = np.round(plot_levels*plot_norm,1).astype(str).tolist()
            cb_labels[-1] = ''  # Remove the final label
            cb.set_ticks(ticks=plot_levels*plot_norm, labels=cb_labels)
        cb.set_label(f'Rainfall ({plot_units})')  # Label the colorbar
        ax.set_title(f"Ensemble standard deviation",size=14)  # This plot's title

        fig.suptitle(f"Jurre Brishti cGAN forecast: Valid {valid_time} - {valid_time + valid_time_delta} UTC")  # Overall title
        plt.tight_layout()  # Looks nicer

        # Save the plot
        if (file_name != None):
            if file_name[-4:] in ['.png','.jpg','.pdf']:
                # If we are making more than one plot
                if (valid_time_start_hour == 'all'):
                    # Append the hour to the file name
                    save_file_name = f"{file_name[:-4]}_{valid_time.hour:02d}{file_name[-4:]}"
                else:  # We are making only one plot
                    save_file_name = file_name  # Use the exact file name specified
                plt.savefig(save_file_name, format=file_name[-3:], bbox_inches='tight')
            else:
                print("ERROR: File type must be specified by '.png', '.jpg' or '.pdf'")

        plt.show()  # Finally draw the plot
     
        
# Plot all ensemble members in the cGAN forecast data at a specified valid time.
# Arguments
#   data                    - An xarray DataSet containing the cGAN rainfall forecasts.
#   valid_time_start_hour=6 - The hour the valid time starts at. Can either be 6, 12, 18 or 0 UTC.
#   style=None              - Options: 'ICPAC', 'ICPAC_heavy', 'KMD', 'EMI', 'EMI_heavy'
#   plot_units='mm/h'       - Can be 'mm/h' (default), 'mm/6h', 'mm/day' or 'mm/week'
#   region='ICPAC'          - can be 'ICPAC', 'Kenya', 'South Sudan', 'Rwanda', 'Burundi', 'Djibouti',
#                             'Eritrea', 'Ethiopia', 'Sudan', 'Somalia', 'Tanzania', 'Uganda'
#   max_num_plots=50        - Maximum number of ensemble members to plot.
#   file_name=None          - If a file name, ending in '.png', '.jpg' or '.pdf' is specified, the
#                             plot is saved in that format.
def plot_GAN_ensemble(data, valid_time_start_hour=6, style=None, plot_units='mm/h', region='ICPAC',
                      max_num_plots=50, file_name=None):
    
    # Get the units to use for plotting
    plot_norm, plot_units = get_plot_normalisation(plot_units)

    # To be consistent with the Harris et. al paper.
    value_range_precip = (0.1, 15 * plot_norm)
    
    # Use a style other than the default
    if (style is not None):
        plot_levels, plot_colours = get_contour_levels(style)

    # Load the border shapefile
    reader = shpreader.Reader("show_forecasts/shapes/GHA_shapes/gha.shp")
    borders_feature = ShapelyFeature(reader.geometries(), ccrs.PlateCarree(), facecolor='none')
    
    if (region != 'ICPAC'):
        
        # Load the regions shapefile
        reader = shpreader.Reader(f"show_forecasts/shapes/{region}_shapes/{region}_region.cpg")
        regions_feature = ShapelyFeature(reader.geometries(), ccrs.PlateCarree(), facecolor='none')
        
        # Get the extent of the region that we are looking at
        region_extent = get_region_extent(region, border_size=0.5)
    
    # Change valid_time_start_hour into the valid_time_idx
    if (str(valid_time_start_hour) == '6'):
        valid_time_idx = 0
    elif (str(valid_time_start_hour) == '12'):
        valid_time_idx = 1
    elif (str(valid_time_start_hour) == '18'):
        valid_time_idx = 2
    elif (str(valid_time_start_hour) == '0'):
        valid_time_idx = 3
    else:
        print("ERROR: valid_time_start_hour must be 6, 12, 18 or 0 hours UTC.")
        return
    
    # Convert the forecast valid time to a datetime.datetime format
    valid_time = datetime64_to_datetime(data['fcst_valid_time'][0,valid_time_idx].values)

    # How many plots will we make
    num_plots = np.min([max_num_plots, data['member'].size])
    num_rows = int(np.ceil(num_plots/5))

    # Define the figure and each axes for the rows and columns
    fig, axs = plt.subplots(nrows=num_rows, ncols=5, subplot_kw={'projection': ccrs.PlateCarree()},
                            figsize=(10,2.8*num_rows+0.9), layout="constrained")

    # axs is a 2 dimensional array of `GeoAxes`. Flatten it into a 1-D array
    axs=axs.flatten()

    # Don't show axes without plots
    for ax_idx in range(num_plots, num_rows*5):
        axs[ax_idx].set_axis_off()

    # For each ensemble member
    for ax_idx in range(num_plots):

        ax=axs[ax_idx]  # First plot (left)
        ax.add_feature(cfeature.COASTLINE, linewidth=1)  # Draw some features to see where we are
        if (region != 'ICPAC'):
            if (region != 'Uganda'):  # Uganda counties are too complicated
                ax.add_feature(regions_feature, linestyle=':')
            ax.set_extent(region_extent, crs=ccrs.PlateCarree())
        ax.add_feature(borders_feature)  # The borders
        ax.add_feature(cfeature.LAKES, linewidth=1,linestyle='-',edgecolor='dimgrey',facecolor='none')
        # Actually make the plot
        if (style == None):
            c = ax.pcolormesh(data['longitude'], data['latitude'], data['precipitation'][0,ax_idx,valid_time_idx,:,:] * plot_norm,
                              norm=colors.LogNorm(*value_range_precip), cmap='YlGnBu', transform=ccrs.PlateCarree())
        else:
            c = ax.contourf(data['longitude'], data['latitude'], data['precipitation'][0,ax_idx,valid_time_idx,:,:] * plot_norm,
                            colors=plot_colours, levels=plot_levels*plot_norm, transform=ccrs.PlateCarree())
        ax.set_title(f"{ax_idx+1}",size=14)  # This plot's title

    # Add a final colorbar with a nice size
    cb = fig.colorbar(c, ax=axs, location='bottom', shrink=0.6, pad=0.01) 
    if (style != None):
        cb_labels = np.round(plot_levels*plot_norm,1).astype(str).tolist()
        cb_labels[-1] = ''  # Remove the final label
        cb.set_ticks(ticks=plot_levels*plot_norm, labels=cb_labels)
    cb.set_label(f'Rainfall ({plot_units})')  # Label the colorbar

    fig.suptitle(f"Jurre Brishti cGAN ensemble: Valid {valid_time} - {valid_time + timedelta(hours=6)} UTC")  # Overall title

    # Save the plot
    if (file_name != None):
        if file_name[-4:] in ['.png','.jpg','.pdf']:
            plt.savefig(file_name, format=file_name[-3:], bbox_inches='tight')
        else:
            print("ERROR: File type must be specified by '.png', '.jpg' or '.pdf'")
    
    plt.show()  # Finally draw the plot


# Plot the chance of rainfall at a rate above a specified threshold.
# Arguments
#   data                    - An xarray DataSet containing the cGAN rainfall forecasts.
#   threshold=2             - We'll plot the chance of rainfall above this threshold rate. The
#                             default is 2 mm/h. The units of threshold is set by plot_units.
#   plot_units='mm/h'       - Can be 'mm/h' (default), 'mm/6h', 'mm/day' or 'mm/week'
#   valid_time_start_hour=6 - The hour the valid time starts at. Can either be 6, 12, 18 or 0 UTC,
#                             or specify 'all' to make all four plots.
#   style=None              - Options: 'ICPAC', 'KMD', 'EMI'
#   show_percentages=False  - Either shows a description (False) or the percentage (True) of
#                             the chance of exceeding the threshold.
#   region='ICPAC'          - Can be 'ICPAC', 'Kenya', 'South Sudan', 'Rwanda', 'Burundi', 'Djibouti',
#                             'Eritrea', 'Ethiopia', 'Sudan', 'Somalia', 'Tanzania', 'Uganda'
#   file_name=None          - If a file name, ending in '.png', '.jpg' or '.pdf' is specified, the
#                             plot is saved in that format.
def plot_GAN_threshold_chance(data, threshold=2, plot_units='mm/h', valid_time_start_hour=6, style=None,
                              show_percentages=False, region='ICPAC', file_name=None):

    # Get the units to use for plotting
    plot_norm, plot_units = get_plot_normalisation(plot_units)

    # Normalise the threshold to the chosen units
    threshold /= plot_norm

    # Plot countour levels and their labels
    plot_levels = [0,2.5,25,50,75,100]  # Percentage chance of rain above this value
    plot_level_names = ['None','Possible','Low','Medium','High','Certain']
    plot_level_percentages = []
    for i in range(len(plot_levels)):
        plot_level_percentages.append(f"{plot_levels[i]}%")

    # Load the border shapefile
    reader = shpreader.Reader("show_forecasts/shapes/GHA_shapes/gha.shp")
    borders_feature = ShapelyFeature(reader.geometries(), ccrs.PlateCarree(), facecolor='none')

    if (region != 'ICPAC'):
                    
        # Load the regions shapefile
        reader = shpreader.Reader(f"show_forecasts/shapes/{region}_shapes/{region}_region.cpg")
        regions_feature = ShapelyFeature(reader.geometries(), ccrs.PlateCarree(), facecolor='none')
        
        # Get the extent of the region that we are looking at
        region_extent = get_region_extent(region, border_size=0.5)

    # Change valid_time_start_hour into the valid_time_idx
    if (str(valid_time_start_hour) == '6'):
        valid_time_idx_list = [0]
    elif (str(valid_time_start_hour) == '12'):
        valid_time_idx_list = [1]
    elif (str(valid_time_start_hour) == '18'):
        valid_time_idx_list = [2]
    elif (str(valid_time_start_hour) == '0'):
        valid_time_idx_list = [3]
    elif (valid_time_start_hour == 'all'):
        valid_time_idx_list = np.arange(len(data['valid_time']))
    else:
        print("ERROR: valid_time_start_hour must be 6, 12, 18 or 0 hours UTC or 'all'.")
        return

    if (len(valid_time_idx_list) == 1):

        # Define the figure and axes
        fig, axs = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(5,5))

        # axs is a `GeoAxes`. Make it into a 1-D array
        axs=[axs]

    if (len(valid_time_idx_list) == 4):
    
        # Define the figure and axes
        fig, axs = plt.subplots(nrows=2, ncols=2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(8,8))

        # axs is a 2 dimensional array of `GeoAxes`. Flatten it into a 1-D array
        axs=axs.flatten()

    # Define the plot colours
    plot_colours = get_threshold_plot_colours(style)

    # There are plots for each valid time
    for idx, valid_time_idx in enumerate(valid_time_idx_list):

        # Convert the forecast valid time to a datetime.datetime format
        valid_time = datetime64_to_datetime(data['fcst_valid_time'][0,valid_time_idx].values)

        # Keep the first valid time for the plot title
        if (idx == 0):
            first_valid_time = valid_time

        # The percentage of ensemble members that exceed the threshold
        plot_data = np.sum(data["precipitation"][0,:,valid_time_idx,:,:] > threshold, axis=0) * 100 / len(data["member"])

        ax = axs[idx]
        ax.gridlines()
        ax.set_facecolor('white')  # For consistency with Harris et. al 2022
        ax.add_feature(cfeature.COASTLINE, linewidth=1)
        if (region != 'ICPAC'):
            if (region != 'Uganda'):  # Uganda counties are too complicated
                ax.add_feature(regions_feature, linestyle=':')
            ax.set_extent(region_extent, crs=ccrs.PlateCarree())
        ax.add_feature(borders_feature)  # The borders
        ax.add_feature(cfeature.LAKES, linewidth=1,linestyle='-',edgecolor='dimgrey',facecolor='none')

        c = ax.contourf(data['longitude'], data['latitude'], plot_data,
                        levels=plot_levels, transform=ccrs.PlateCarree(), colors=plot_colours)
        ax.set_title(f"Valid {valid_time.time()} - {(valid_time + timedelta(hours=6)).time()}", size=14)

        cb = plt.colorbar(c, fraction=0.04)
        #cb.ax.tick_params(labelsize=18)
        if show_percentages:
            #cb.set_label(f'% chance',size=18)
            cb.set_ticks(ticks=plot_levels,labels=plot_level_percentages)
        else:
            cb.set_ticks(ticks=plot_levels,labels=plot_level_names)

    title_string = f"""Jurre Brishti cGAN ensemble: {first_valid_time.date()} - {(valid_time + timedelta(hours=6)).date()} UTC
    Chance of rainfall above {threshold*plot_norm:.1f} {plot_units}."""
    fig.suptitle(title_string)  # Overall title
    plt.tight_layout()  # Looks nicer

    # Save the plot
    if (file_name != None):
        if file_name[-4:] in ['.png','.jpg','.pdf']:
            plt.savefig(file_name, format=file_name[-3:], bbox_inches='tight')
        else:
            print("ERROR: File type must be specified by '.png', '.jpg' or '.pdf'")

    plt.show()


# Draw a marker at a specific location.
#   name=None               - An optional string containing the name of the point.
#   latitude=None           - The latitude of the point.
#   longitude=None          - The longitude of the point.
#   region='ICPAC'          - Can be 'ICPAC', 'Kenya', 'South Sudan', 'Rwanda', 'Burundi', 'Djibouti',
#                             'Eritrea', 'Ethiopia', 'Sudan', 'Somalia', 'Tanzania', 'Uganda'
#   file_name=None          - If a file name, ending in '.png', '.jpg' or '.pdf' is specified, the
#                             plot is saved in that format.
def plot_location_marker(location_name=None, latitude=None, longitude=None, region='ICPAC', file_name=None):

    if (((latitude == None) and (longitude != None)) or
        ((latitude != None) and (longitude == None))):
        print("ERROR: Either don't specify latitude and longitude or specify both.")
        return

    # lattiude and longitude are not specified
    if (latitude == None) and (longitude == None):
        # Select the location from the location name
        location_found = False
        for location in locations:
            if ((location["name"] == location_name) and ((location["country"] == region) or (region == 'ICPAC'))):
                location_found = True
                break
        if (not location_found):
            print(f"ERROR: Location '{location_name}' is not in the list of locations.")
            if (region != 'ICPAC'):
                # A country is specified
                print(f"       Searching in {region}. Perhaps try a different region.")
            return
    
    else:  # latitude and longitude are specified
        location = {"name":location_name,
                    "country":"",
                    "latitude":latitude,
                    "longitude":longitude}

    # Load the border shapefile
    reader = shpreader.Reader("show_forecasts/shapes/GHA_shapes/gha.shp")
    borders_feature = ShapelyFeature(reader.geometries(), ccrs.PlateCarree(), facecolor='none')

    if (region != 'ICPAC') and (region != None):
                    
        # Load the regions shapefile
        reader = shpreader.Reader(f"show_forecasts/shapes/{region}_shapes/{region}_region.cpg")
        regions_feature = ShapelyFeature(reader.geometries(), ccrs.PlateCarree(), facecolor='none')
        
    # Get the extent of the region that we are looking at
    region_extent = get_region_extent(region, border_size=0.5)

    # Check that the point specified is within the region specified
    if (pt_in_rect([location['longitude'], location['latitude']], region_extent) != True):
        print(f'ERROR: Longitude and latitude specified is not in {region}.')
        return

    # Define the figure and axes
    fig, ax = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(5,5))

    ax.gridlines()
    ax.set_facecolor('white')  # For consistency with Harris et. al 2022
    ax.add_feature(cfeature.COASTLINE, linewidth=1)
    if (region != 'ICPAC'):
        if (region != 'Uganda'):  # Uganda counties are too complicated
            ax.add_feature(regions_feature, linestyle=':')
    ax.set_extent(region_extent, crs=ccrs.PlateCarree())
    ax.add_feature(borders_feature)  # The borders
    ax.add_feature(cfeature.LAKES, linewidth=1,linestyle='-',edgecolor='dimgrey',facecolor='none')

    # Plot the location
    plt.plot(location['longitude'], location['latitude'], 'ro')

    if (location_name != None):
        if (location["country"] != ""):
            ax.set_title(f"{location_name}, {location["country"]}", size=14)
        else:
            ax.set_title(f"{location_name}", size=14)

    plt.tight_layout()  # Looks nicer

    # Save the plot
    if (file_name != None):
        if file_name[-4:] in ['.png','.jpg','.pdf']:
            plt.savefig(file_name, format=file_name[-3:], bbox_inches='tight')
        else:
            print("ERROR: File type must be specified by '.png', '.jpg' or '.pdf'")

    plt.show()


# Plot histograms from the ensemble values of rainfall at a specified location.
# Arguments
#   location_name     - The location name corresponding to a name in the list returned by print_locations().
#   country=None      - For the case of the same name in different countries, the country can optionally
#                       be selected from 'Kenya', 'South Sudan', 'Rwanda', 'Burundi', 'Djibouti',
#                       'Eritrea', 'Ethiopia', 'Sudan', 'Somalia', 'Tanzania', 'Uganda'.
#   latitude=None     - Optional latiude in degrees north. If specified along with a longitude, this
#                       location is used instead of any named location in the list returned by
#                       print_locations().
#   longitude=None    - Optional longitude in degrees north. If specified along with a latitude, this
#                       location is used instead of any named location in the list returned by
#                       print_locations().
#   plot_units='mm/h' - Can be 'mm/h' (default), 'mm/6h', 'mm/day' or 'mm/week'
#   num_bins=10       - Number of evenly spaced bins in the histogram
#   probability=None  - Plot a line indicating the amount of rain that will be exceeded with a given probability.
def plot_GAN_local_histograms(data, location_name, country=None, latitude=None, longitude=None, plot_units='mm/h',
                              num_bins=10, probability=None):

    if (((latitude == None) and (longitude != None)) or
        ((latitude != None) and (longitude == None))):
        print("ERROR: Either don't specify latitude and longitude or specify both.")
        return

    # lattiude and longitude are not specified
    if (latitude == None) and (longitude == None):
        # Select the location from the location name
        location_found = False
        for location in locations:
            if ((location["name"] == location_name) and ((location["country"] == country) or (country == None))):
                location_found = True
                break
        if (not location_found):
            print(f"ERROR: Location '{location_name}' is not in the list of locations.")
            return
    
    else:  # latitude and longitude are specified
        location = {"name":location_name,
                    "country":"",
                    "latitude":latitude,
                    "longitude":longitude}

    # Select the nearest longitude and latitude indices
    lat_idx = np.argmin(np.abs(location["latitude"] - data["latitude"].values))
    lon_idx = np.argmin(np.abs(location["longitude"] - data["longitude"].values))

    # Get the units to use for plotting
    plot_norm, plot_units = get_plot_normalisation(plot_units)

    # Get the avaliable valid times
    valid_time_idx_list = np.arange(len(data['valid_time']))

    # Define the figure and axes
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8,8))

    # axs is a 2 dimensional array of `GeoAxes`. Flatten it into a 1-D array
    axs=axs.flatten()

    # There are plots for each valid time
    for idx, valid_time_idx in enumerate(valid_time_idx_list):

        # Convert the forecast valid time to a datetime.datetime format
        valid_time = datetime64_to_datetime(data['fcst_valid_time'][0,valid_time_idx].values)

        # Keep the first valid time for the plot title
        if (idx == 0):
            first_valid_time = valid_time

        # The data to plot
        plot_data = data['precipitation'][0,:,valid_time_idx,lat_idx,lon_idx].values * plot_norm

        # Plot the histogram
        ax = axs[idx]
        style = {'facecolor': 'tab:blue', 'edgecolor': 'black', 'linewidth': 1}
        ax.hist(plot_data, num_bins, **style)

        # The percentile specified
        if (probability != None):
            # Convert the probability of excedence to a percentile and the corresponding percentile index
            data_sorted = np.sort(plot_data)
            percentile = 1-probability
            percentile_idx = int(np.floor(len(data_sorted) * percentile))
            exceedance_value = np.sort(plot_data)[percentile_idx]
            ax.axvline(x=exceedance_value, color='r', label=f"{percentile*100:0g}%: {exceedance_value:.3g} {plot_units}")
            ax.legend()

        ax.grid()
        ax.set_xlabel(f'Rainfall ({plot_units})')
        ax.set_ylabel('Number of ensemble members')
        ax.set_title(f"Valid {valid_time.time()} - {(valid_time + timedelta(hours=6)).time()} UTC")

    title_string = f"""Jurre Brishti cGAN ensemble {first_valid_time.date()} - {(valid_time + timedelta(hours=6)).date()} UTC
    {location["name"]}, {location["country"]} ({location["latitude"]:.4f}N, {location["longitude"]:.4f}E)"""
    plt.suptitle(title_string)  # Overall title    
    plt.tight_layout()  # Looks nicer
    plt.show()


# Plot the amount of rain that will be exceeded with a given probability.
# Arguments
#   data_sorted             - An xattay Dataset containing the cGAN rainfall forecasts sorted along the
#                             ensemble axis.
#   probability=0.05        - We'll plot the amount of rainfall at this threshold probability. The
#                             default is 0.05 or 1/20. So by default there is a 1 in 20 chance of
#                             rainfall above this value.
#   valid_time_start_hour=6 - The hour the valid time starts at. Can either be 6, 12, 18 or 0 UTC,
#                             or specify 'all' to make all four plots.
#   style=None              - Options: 'ICPAC', 'KMD', 'EMI'
#   plot_units='mm/h'       - Can be 'mm/h' (default), 'mm/6h', 'mm/day' or 'mm/week'
#   region='ICPAC'          - Can be 'ICPAC', 'Kenya', 'South Sudan', 'Rwanda', 'Burundi', 'Djibouti',
#                             'Eritrea', 'Ethiopia', 'Sudan', 'Somalia', 'Tanzania', 'Uganda'
#   file_name=None          - If a file name, ending in '.png', '.jpg' or '.pdf' is specified, the
#                             plot is saved in that format.
def plot_GAN_exceedance(data_sorted, probability=0.05, valid_time_start_hour=6, style=None,
                        plot_units='mm/h', region='ICPAC', file_name=None):

    # Get the units to use for plotting
    plot_norm, plot_units = get_plot_normalisation(plot_units)
    
    # To be consistent with the Harris et. al paper.
    value_range_precip = (0.1, 15 * plot_norm)
    
    # Use a style other than the default
    if (style is not None):
        plot_levels, plot_colours = get_contour_levels(style)

    # Load the border shapefile
    reader = shpreader.Reader("show_forecasts/shapes/GHA_shapes/gha.shp")
    borders_feature = ShapelyFeature(reader.geometries(), ccrs.PlateCarree(), facecolor='none')

    if (region != 'ICPAC'):
                    
        # Load the regions shapefile
        reader = shpreader.Reader(f"show_forecasts/shapes/{region}_shapes/{region}_region.cpg")
        regions_feature = ShapelyFeature(reader.geometries(), ccrs.PlateCarree(), facecolor='none')
        
        # Get the extent of the region that we are looking at
        region_extent = get_region_extent(region, border_size=0.5)

    # Change valid_time_start_hour into the valid_time_idx
    if (str(valid_time_start_hour) == '6'):
        valid_time_idx_list = [0]
    elif (str(valid_time_start_hour) == '12'):
        valid_time_idx_list = [1]
    elif (str(valid_time_start_hour) == '18'):
        valid_time_idx_list = [2]
    elif (str(valid_time_start_hour) == '0'):
        valid_time_idx_list = [3]
    elif (valid_time_start_hour == 'all'):
        valid_time_idx_list = np.arange(len(data_sorted['valid_time']))
    else:
        print("ERROR: valid_time_start_hour must be 6, 12, 18 or 0 hours UTC or 'all'.")
        return

    # Convert the probability of excedence to a percentile and the corresponding percentile index
    percentile = 1-probability
    percentile_idx = int(np.floor(len(data_sorted["member"]) * percentile))

    if (len(valid_time_idx_list) == 1):

        # Define the figure and axes
        fig, axs = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(5,5))

        # axs is a `GeoAxes`. Make it into a 1-D array
        axs=[axs]

    if (len(valid_time_idx_list) == 4):
    
        # Define the figure and axes
        fig, axs = plt.subplots(nrows=2, ncols=2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(8,8))

        # axs is a 2 dimensional array of `GeoAxes`. Flatten it into a 1-D array
        axs=axs.flatten()

    # There are plots for each valid time
    for idx, valid_time_idx in enumerate(valid_time_idx_list):

        # Convert the forecast valid time to a datetime.datetime format
        valid_time = datetime64_to_datetime(data_sorted['fcst_valid_time'][0,valid_time_idx].values)

        # Keep the first valid time for the plot title
        if (idx == 0):
            first_valid_time = valid_time

        # The rainfall at the given probability
        data_to_plot = data_sorted["precipitation"][0,percentile_idx,valid_time_idx,:,:] * plot_norm

        ax = axs[idx]
        ax.gridlines()
        ax.set_facecolor('white')  # For consistency with Harris et. al 2022
        ax.add_feature(cfeature.COASTLINE, linewidth=1)
        if (region != 'ICPAC'):
            if (region != 'Uganda'):  # Uganda counties are too complicated
                ax.add_feature(regions_feature, linestyle=':')
            ax.set_extent(region_extent, crs=ccrs.PlateCarree())
        ax.add_feature(borders_feature)  # The borders
        ax.add_feature(cfeature.LAKES, linewidth=1,linestyle='-',edgecolor='dimgrey',facecolor='none')

        # Actually make the plot
        if (style == None):
            c = ax.pcolormesh(data_sorted['longitude'], data_sorted['latitude'], data_to_plot,
                              norm=colors.LogNorm(*value_range_precip), cmap='YlGnBu', transform=ccrs.PlateCarree())
            cb = plt.colorbar(c, fraction=0.04)  # Add a colorbar with a nice size
        else:
            c = ax.contourf(data_sorted['longitude'], data_sorted['latitude'], data_to_plot,
                            colors=plot_colours, levels=plot_levels*plot_norm, transform=ccrs.PlateCarree())
            cb = plt.colorbar(c, fraction=0.04)  # Add a colorbar with a nice size
            cb_labels = np.round(plot_levels*plot_norm,1).astype(str).tolist()
            cb_labels[-1] = ''  # Remove the final label
            cb.set_ticks(ticks=plot_levels*plot_norm, labels=cb_labels)
        cb.set_label(f'Rainfall ({plot_units})')  # Label the colorbar
        ax.set_title(f"Valid {valid_time.strftime("%H:00")} - {(valid_time + timedelta(hours=6)).strftime("%H:00")} UTC", size=14)

    title_string = f"""Jurre Brishti cGAN ensemble: Valid {first_valid_time.date()} - {(valid_time + timedelta(hours=6)).date()}
    {percentile*100:0g}% chance that rainfall is below this value."""
    fig.suptitle(title_string)  # Overall title
    plt.tight_layout()  # Looks nicer

    # Save the plot
    if (file_name != None):
        if file_name[-4:] in ['.png','.jpg','.pdf']:
            plt.savefig(file_name, format=file_name[-3:], bbox_inches='tight')
        else:
            print("ERROR: File type must be specified by '.png', '.jpg' or '.pdf'")

    plt.show()
