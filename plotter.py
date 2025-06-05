import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re

class Plotter:
    def __init__(self, colors):
        self.colour_palette = colors
    
    def plotBoxplot(self, data_df, measurement_name, y_axis_names_dic):
        data_df = data_df.reset_index(drop=True)
        if data_df.empty:
            print(f"Measurement {measurement_name} not found!")
        else:
            data_df["value"] = pd.to_numeric(data_df["value"], errors='coerce')
            data_df = data_df.dropna(subset = ["value"])
            plt.figure(figsize=(6, 4))
            sns.set_theme(style="whitegrid")
            ax_box = sns.boxplot(data = data_df, hue = "variable", x = "variable", y = "value", palette="magma") 
            ax_box.set_title(measurement_name, fontweight = "bold")
            ax_box.set(xlabel = "", ylabel = y_axis_names_dic[measurement_name])
            ax_box.tick_params(labelrotation = 30)
            sns.stripplot(x = "variable", y = "value", hue = "variable", data=data_df, palette='dark:black', size=3, ax = ax_box)
            plt.xlabel(None)
            plt.show()
            #plt.savefig(f"{measurement_name}", format="svg")

    def plotViolinplot(self, data_df, measurement_name, y_axis_names_dic):
        if data_df.empty:
            print(f"Measurement {measurement_name} not found!")
        else:
            plt.figure(figsize=(6, 4))
            sns.set_theme(style="whitegrid")
            ax_vio = sns.violinplot(data = data_df, x = "variable", y = "value", hue= "variable", palette="magma")
            ax_vio.set_title(measurement_name)
            ax_vio.set(xlabel = "", ylabel = y_axis_names_dic[measurement_name])
            sns.stripplot(x="variable", y = "value", data=data_df, palette='dark:black', size=3, ax=ax_vio)
            plt.xlabel(None)
            plt.show()        
            #plt.savefig(f"{measurement_name}", format="svg")

    def superplot5(self, combined, measurement, y_axis_names_dic, stat_list, directory_for_saving, des_stat_param, save_SVG_param, save_PNG_param):
        print(stat_list)
        if stat_list is not None:
            stat_measurements = [dic for dic in stat_list if dic is not None and dic.get("measurement") == measurement]
        else: 
            stat_measurements = []
             
        sns.set_theme(style="whitegrid")
        combined = combined.reset_index(drop = True)
        combined[measurement] = pd.to_numeric(combined[measurement], errors='coerce')
        combined = combined.dropna(subset=[measurement])

        if(not combined[measurement].empty):
            ReplicateAverages = combined.groupby(['Group', 'Replicate'], as_index=False).agg({str(measurement) : "mean"}).reset_index()

            filename = re.sub(r'[\\/*?:"<>|]', "_", measurement)
            write_path = directory_for_saving + f"/{filename}.csv"
            with open(write_path, "w") as csv: 
                csv.write(f"Plotting results for {measurement}\n")
                csv.write(f"{os.path.abspath(write_path)} \t")
                csv.write("\n")
                combined[[measurement, 'Group', 'Replicate']].to_csv(write_path, mode = "a", sep=",", index = False)
                      
            fig, ax = plt.subplots()
            if measurement in ["orientation vector x [micron] (vector from first to last skeleton point)" , 
                "orientation vector y [micron] (vector from first to last skeleton point)",
                "orientation vector z [micron] (vector from first to last skeleton point)",
                "Maximum span [micron]"]:
                sns.violinplot(x="Group", y=measurement, data=combined, ax=ax, color="lightgrey", inner="quart", alpha = 0.6, density_norm='count')
            else:
                sns.violinplot(x="Group", y=measurement, data=combined, ax=ax, color="lightgrey", inner="quart", alpha = 0.6, density_norm='count', cut = 0)

            #sns.violinplot(x="Group", y=measurement, data=combined, ax=ax, color="lightgrey", inner="quart", alpha = 0.6, density_norm='count')
            sns.stripplot(x="Group", y=measurement, hue="Replicate", data=combined, alpha=0.7, edgecolor="darkgray", linewidth=0.8)
            sns.swarmplot(x = "Group", y=measurement, hue="Replicate", size=15, edgecolor="k", linewidth=2, data=ReplicateAverages, alpha=0.8, ax = ax)
            ax.legend_.remove()
            ax.set_title(measurement, weight = "bold")
            ax.set(xlabel = "", ylabel = y_axis_names_dic[measurement])
            ax.tick_params(labelrotation = 0)
            increment_len = len(des_stat_param[des_stat_param["Measurement"] == measurement]["Group"].unique())

            if stat_measurements:
                for stat_measurement in stat_measurements:
                    if stat_measurement["pairwise"] is not None:
                        label_performed_test = f"Statistical test: {stat_measurement['hypothesis test']}"
                        x_ticks_group = ax.get_xticks()
                        xticklabels = [tick.get_text() for tick in ax.get_xticklabels()]
                        group_positions = {grp: pos for grp, pos in zip(xticklabels, x_ticks_group)}
                        max_y_val = combined.groupby("Group")[measurement].max() 
                        overall_range = combined[measurement].max() - combined[measurement].min()
                        offset = overall_range * 0.10
                        max_y = combined[measurement].max()

                        for increment, (comparison, values) in enumerate (stat_measurement["pairwise"].items()):
                            try:
                                group1, group2 = comparison.split(" vs ")
                                label_performed_test = label_performed_test + f"\n {group1} vs. {group2}: {values[0]} (p-value: {round(values[1], 5)})" 
                            except ValueError:
                                print(f"Comparison '{comparison}' not formatted as 'Group1 vs Group2'. Skipping this comparison.")
                                continue  
                            except KeyError:
                                print(f"Comparison '{comparison}' not formatted as 'Group1 vs Group2'. Skipping this comparison.")

                            x1 = group_positions.get(group1)
                            x2 = group_positions.get(group2)
                            if x1 is None or x2 is None:
                                print(f"Could not find x positions for comparison '{comparison}'.")

                            y1 = max_y_val.get(group1, 0)
                            y2 = max_y_val.get(group2, 0)
                            y = max_y + offset * (1.0 + 2.0*increment)

                            #add padding above statistics bar
                            current_ax_lim = ax.get_ylim()  # returns (lower, upper)
                            lower, upper = current_ax_lim
                            new_ax_lim = upper + offset * 2  # offset is already a float
                            ax.set_ylim(lower, new_ax_lim)
                            ax.plot([x1, x2], [y, y], lw=1.5, c='black')
                            ax.plot([x1, x1], [y, y - offset/2], lw=1.5, c='black')
                            ax.plot([x2, x2], [y, y - offset/2], lw=1.5, c='black')
                            ax.text((x1 + x2) / 2, y + offset/10, values[0], ha='center', va='bottom', color='black', fontsize=15)
            else:
                label_performed_test = f"Statistical test: No statistical test was performed."
                #fig.suptitle(label_performed_test, x=0.01, y=-0.03, ha='left', va='bottom', fontsize=10)
            label_performed_test = label_performed_test + "\n"
            for i, row in des_stat_param[des_stat_param["Measurement"] == measurement].iterrows():
                label_performed_test = label_performed_test + f"\n {row.loc['Group']}: Mean = {row.loc['Mean']} (±{row.loc['Std']}), count = {row.loc['Count']}"

            subtitle = fig.suptitle(label_performed_test, x=0.01, y=-0.1 * increment_len, ha='left', va='bottom', fontsize=10)
            
            if save_PNG_param == True: 
                plt.savefig(f"{directory_for_saving}/{filename}.png", dpi = 300, bbox_inches='tight', bbox_extra_artists = [subtitle])
                print(f"now, I would save the figure under  {directory_for_saving}/{filename}.png")
            if save_SVG_param == True: 
                plt.savefig(f"{directory_for_saving}/{filename}.svg")
                print(f"now, I would save the figure under  {directory_for_saving}/{filename}.svg")
            plt.show()
            plt.close(fig)
        else:
            print(measurement + " was not found")
