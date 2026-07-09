    bin_labels = [np.mean([dist_bins[i], dist_bins[i+1]]) for i in range(len(dist_bins)-1)]
    dict_
    for _, dict_cfg in cfg.items():
        df_crstar = df['cruise'] == dict_cfg['cruise']


    # 3. Bin Latitudes
    # Use the midpoints of the bins to plot against, or use bin strings. 
    # For a continuous line plot, midpoints work best.
    df['lat_bin_center'] = pd.cut(df['latitude'], bins=lat_bins, labels=bin_labels, include_lowest=True)
    
    # Drop rows that fall outside the specified latitude bins
    df = df.dropna(subset=['lat_bin_center'])
    df['lat_bin_center'] = df['lat_bin_center'].astype(float)

    # 4. Aggregate Data (Calculate Mean and Standard Deviation)
    # Grouping by the bin center, cruise, and target
    agg_df = df.groupby(['lat_bin_center', 'cruise', 'target']).agg(
        mean_concentration=('concentration', 'mean'),
        std_concentration=('concentration', 'std'),
        count=('concentration', 'count')
    ).reset_index()

    # Sort by latitude bin center to ensure lines connect chronologically/geographically
    agg_df = agg_df.sort_values(by='lat_bin_center')

    # 5. Plotting
    plt.figure(figsize=(10, 6))
    
    # Get unique cruises to loop through for custom error bar rendering
    cruises = agg_df['cruise'].unique()
    # Use a standard color palette for consistency across lines and error bars
    colors = sns.color_palette("Set1", len(cruises))
    
    for cruise, color in zip(cruises, colors):
        cruise_data = agg_df[agg_df['cruise'] == cruise]
        
        # Fill NaN std deviations with 0 (happens if a bin has only 1 data point)
        y_err = cruise_data['std_concentration'].fillna(0).values
        
        plt.errorbar(
            x=cruise_data['lat_bin_center'],
            y=cruise_data['mean_concentration'],
            yerr=y_err,
            label=f"Cruise: {cruise}",
            fmt='-o', # line (-) and dot (o)
            color=color,
            capsize=4,
            markersize=6,
            linewidth=2
        )

    plt.xlabel('Latitude (Bin Centers)')
    plt.ylabel(f'Mean Concentration of {target_col}')
    plt.title(f'{target_col} Concentration Profile by Cruise and Latitude')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()
    
    # Save the plot
    output_plot_path = f"{plot_base}.png"
    plt.tight_layout()
    plt.savefig(output_plot_path, dpi=300)
    click.echo(f"Successfully saved plot to: {output_plot_path}")
